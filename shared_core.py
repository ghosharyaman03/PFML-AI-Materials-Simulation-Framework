#!/usr/bin/env python3
"""
shared_core.py
===============
Common code shared by all three pipeline stages (translator, builder,
validator) plus the orchestrator. Nothing in here talks to the network,
Ollama, or pycalphad at import time -- those stay lazy so each stage can
import this module cheaply.

Contains, unchanged from the original monolithic run_agents.py:
  - EquilibriumSolver          (dynamic thermodynamic value computation)
  - Configuration schema       (SECTIONS, PARAM_ALIASES, MULTI_LINE_KEYS)
  - Template engine            (generate/parse .in files)
  - discover_systems           (scans Example_Systems/ for presets)
  - Display helpers            (used by the orchestrator's console output)
"""

import os
import re
import json
import copy
import math
import glob

# ============================================================================
# THERMODYNAMIC SOLVER (pycalphad-based)
# ============================================================================

class EquilibriumSolver:
    """Compute ceq/cfill/c_guess/slopes via pycalphad for N-phase simulations.

    Only activated when the requested phase count exceeds the pre-defined
    template data.  All calculations happen in Python before the C binary
    is invoked — the binary itself is never modified.
    """

    _ELEMENT_REFS = {
        'AL': 'FCC_A1', 'ZN': 'HCP_ZN', 'NI': 'FCC_A1', 'CU': 'FCC_A1',
        'MO': 'BCC_A2', 'CR': 'BCC_A2', 'FE': 'BCC_A2', 'TI': 'HCP_A3',
        'NB': 'BCC_A2', 'CO': 'HCP_A3', 'W':  'BCC_A2', 'TA': 'BCC_A2',
    }

    def __init__(self, base_dir="/app/microsim_gp"):
        self._base_dir = base_dir

    # ------------------------------------------------------------------
    # TDB management
    # ------------------------------------------------------------------

    def find_tdb(self, tdb_name):
        """Search known directories for *tdb_name*.  Returns path or None."""
        if not tdb_name:
            return None
        for sub in ("", "tdbs"):
            candidate = (os.path.join(self._base_dir, sub, tdb_name) if sub
                         else os.path.join(self._base_dir, tdb_name))
            if os.path.isfile(candidate):
                return candidate
        return None

    def find_best_tdb(self, components, needed_phases):
        """Scan every .tdb file under <base_dir>/tdbs and return the path
        of the first one that defines both the required elements AND the
        required phase names (e.g. L12_FCC). Returns None if nothing in
        the folder covers this system.
        """
        from pycalphad import Database
        tdb_dir = os.path.join(self._base_dir, "tdbs")
        if not os.path.isdir(tdb_dir):
            return None
        needed_elements = {c.upper() for c in components}
        needed_phase_names = {p.upper() for p in needed_phases}
        for fname in sorted(os.listdir(tdb_dir)):
            if not fname.lower().endswith(".tdb"):
                continue
            candidate_path = os.path.join(tdb_dir, fname)
            try:
                db = Database(candidate_path)
                db_elements = {str(e).upper() for e in db.elements}
                db_phases = {str(p).upper() for p in db.phases.keys()}
                if needed_elements.issubset(db_elements) and needed_phase_names.issubset(db_phases):
                    print(f"  [EquilibriumSolver] Found matching TDB: {candidate_path}")
                    return candidate_path
            except Exception:
                continue
        return None

    def generate_base_tdb(self, components, tdb_name="auto_generated.tdb"):
        """Create a minimal TDB so pycalphad has *something* to work with.

        The generated database is intentionally simple (ideal-solution
        with a single interaction parameter).  Replace it with a properly
        assessed TDB for production simulations.
        """
        tdb_dir = os.path.join(self._base_dir, "tdbs")
        os.makedirs(tdb_dir, exist_ok=True)
        tdb_path = os.path.join(tdb_dir, tdb_name)

        lines = [
            "$ Auto-generated base TDB for EquilibriumSolver",
            "$ Replace with a properly assessed TDB for production use",
            "",
        ]
        for c in components:
            ref = self._ELEMENT_REFS.get(c.upper(), 'FCC_A1')
            lines.append(f"ELEMENT {c.upper()} {ref} 0 0 0 !")
        lines.append("ELEMENT VA VACUUM 0 0 0 !")
        lines.append("")

        comp_str = ",".join(c.upper() for c in components)
        for phase_name in ("LIQUID", "FCC_A1", "BCC_A2"):
            lines.append(f"PHASE {phase_name} % 1 1 !")
            lines.append(f"CONSTITUENT {phase_name} : {comp_str} : !")
            for c in components:
                lines.append(f"PARAMETER G({phase_name},{c.upper()};0) 0 0; 6000 N !")
            for i in range(len(components)):
                for j in range(i + 1, len(components)):
                    c1, c2 = components[i].upper(), components[j].upper()
                    lines.append(
                        f"PARAMETER G({phase_name},{c1},{c2};0) 0 -10000; 6000 N !"
                    )
            lines.append("")

        with open(tdb_path, "w") as f:
            f.write("\n".join(lines) + "\n")
        print(f"  [EquilibriumSolver] Generated base TDB: {tdb_path}")
        return tdb_path

    # ------------------------------------------------------------------
    # Equilibrium calculation
    # ------------------------------------------------------------------

    @staticmethod
    def _extract_phase_compositions(eq_result, phase_name, indep_comps):
        """Pull median equilibrium mole-fractions from a pycalphad result."""
        import numpy as np
        comp_values = []
        for comp in indep_comps:
            x_arr = (eq_result.X.sel(component=comp)
                     .where(eq_result.Phase == phase_name)
                     .values.flatten())
            x_valid = x_arr[
                np.isfinite(x_arr) & (x_arr > 1e-12) & (x_arr < 1.0 - 1e-12)
            ]
            if len(x_valid) > 0:
                comp_values.append(round(float(np.median(x_valid)), 10))
            else:
                comp_values.append(None)
        return comp_values

    def compute_params(self, config, n_phases):
        """Compute ceq, cfill, c_guess (and optionally slopes) for *n_phases*.

        Parameters
        ----------
        config : dict
            The current simulation config (must already contain COMPONENTS,
            NUMCOMPONENTS, T, tdbfname, phase_map, etc.).
        n_phases : int
            Target number of phases.

        Returns
        -------
        dict  with keys ``'ceq'``, ``'cfill'``, ``'c_guess'`` — each a list
        of ``[phase_a, phase_b, val1, ...]`` rows ready to be stored in
        *config*.
        """
        try:
            import numpy as np
            from pycalphad import Database, equilibrium
            import pycalphad.variables as v
        except ImportError:
            raise ImportError(
                "pycalphad is required for dynamic phase computation. "
                "Install with: pip install pycalphad numpy"
            )

        n_comp = int(config.get('NUMCOMPONENTS', 2))
        temp   = float(config.get('T', 1000))
        tdb_name = str(config.get('tdbfname', ''))

        # --- Components --------------------------------------------------
        components_raw = config.get('COMPONENTS', [])
        if isinstance(components_raw, str):
            components_raw = [components_raw]
        components = [c.strip().upper() for c in components_raw]

        # --- Phase map (sim-phase -> TDB-phase) --------------------------
        phase_map_raw = config.get('phase_map', [])
        if isinstance(phase_map_raw, str):
            phase_map_raw = [phase_map_raw]
        phase_map = [p.strip() for p in phase_map_raw]
        while len(phase_map) < n_phases:
            phase_map.append(phase_map[-1] if phase_map else 'FCC_A1')

        # --- TDB ----------------------------------------------------------
        unique_tdb_phases = list(set(phase_map))   # phases this run actually needs

        tdb_path = self.find_tdb(tdb_name)
        if tdb_path:
            # tdbfname in preset Input.in files is often stale/copy-pasted
            # metadata -- verify the database actually defines the elements
            # AND phases this system needs (not just the elements).
            try:
                candidate_db = Database(tdb_path)
                db_elements = {str(e).upper() for e in candidate_db.elements}
                db_phases = {str(p).upper() for p in candidate_db.phases.keys()}
                needed_elements = {c.upper() for c in components}
                needed_phase_names = {p.upper() for p in unique_tdb_phases}
                missing_elements = needed_elements - db_elements
                missing_phases = needed_phase_names - db_phases
                if missing_elements or missing_phases:
                    reason = []
                    if missing_elements:
                        reason.append(f"elements {sorted(missing_elements)}")
                    if missing_phases:
                        reason.append(f"phases {sorted(missing_phases)}")
                    print(f"  [EquilibriumSolver] '{tdb_name}' is missing "
                          f"{' and '.join(reason)} -- this preset's tdbfname "
                          f"field is stale/incorrect. Ignoring it.")
                    tdb_path = None
            except Exception as e:
                print(f"  [EquilibriumSolver] Could not read '{tdb_name}' ({e}). "
                      f"Ignoring it.")
                tdb_path = None

        if not tdb_path:
            # Search every .tdb already sitting in microsim_gp/tdbs/ for one
            # that actually covers the elements AND phases this run needs,
            # instead of jumping straight to the bare-bones fallback.
            tdb_path = self.find_best_tdb(components, unique_tdb_phases)

        if not tdb_path:
            print(f"  [EquilibriumSolver] No TDB in tdbs/ covers "
                  f"{components} with phases {unique_tdb_phases}. "
                  f"Generating a minimal fallback (LIQUID/FCC_A1/BCC_A2 only).")
            tdb_path = self.generate_base_tdb(
                components, tdb_name or "auto_generated.tdb"
            )
        print(f"  [EquilibriumSolver] Using TDB : {tdb_path}")
        print(f"  [EquilibriumSolver] Computing equilibrium at T={temp} K "
              f"for {n_phases} phases")

        db = Database(tdb_path)
        pycalphad_comps = list(set(components + ['VA']))
        indep_comps = components[:-1]            # last component is dependent

        # --- Equilibrium for each unique TDB-phase pair -------------------
        unique_tdb_phases = list(set(phase_map))
        eq_compositions = {}                     # (tp_a,tp_b) -> {tp: [x1,...]}
        failed_pairs = []

        for tp_a in unique_tdb_phases:
            for tp_b in unique_tdb_phases:
                pair_key = (tp_a, tp_b)
                if pair_key in eq_compositions:
                    continue

                calc_phases = list(set([tp_a, tp_b]))
                try:
                    if n_comp == 2:
                        conds = {
                            v.T: temp, v.P: 101325, v.N: 1,
                            v.X(indep_comps[0]): np.linspace(0.005, 0.995, 80),
                        }
                    elif n_comp == 3:
                        conds = {
                            v.T: temp, v.P: 101325, v.N: 1,
                            v.X(indep_comps[0]): np.linspace(0.005, 0.95, 40),
                            v.X(indep_comps[1]): np.linspace(0.005, 0.45, 20),
                        }
                    else:
                        conds = {v.T: temp, v.P: 101325, v.N: 1}
                        for ic in indep_comps:
                            conds[v.X(ic)] = np.linspace(0.01, 0.5, 15)

                    eq_result = equilibrium(db, pycalphad_comps,
                                            calc_phases, conds)

                    pair_data = {}
                    for pname in calc_phases:
                        pair_data[pname] = self._extract_phase_compositions(
                            eq_result, pname, indep_comps
                        )
                    eq_compositions[pair_key] = pair_data
                    if tp_a != tp_b:
                        eq_compositions[(tp_b, tp_a)] = pair_data

                except Exception as e:
                    print(f"  [EquilibriumSolver] Warning: calc failed for "
                          f"({tp_a}, {tp_b}): {e}")
                    failed_pairs.append((tp_a, tp_b))
                    fallback = [1.0 / n_comp] * (n_comp - 1)
                    pair_data = {p: list(fallback) for p in calc_phases}
                    eq_compositions[pair_key] = pair_data
                    if tp_a != tp_b:
                        eq_compositions[(tp_b, tp_a)] = pair_data

        # --- Build ceq / cfill / c_guess rows ----------------------------
        ceq_entries   = []
        cfill_entries = []
        cguess_entries = []
        fallback_comp = [1.0 / n_comp] * (n_comp - 1)

        for i in range(n_phases):
            for j in range(n_phases):
                tp_i = phase_map[i]
                tp_j = phase_map[j]

                # NOTE: previously this skipped (i, j) pairs when two
                # different sim-phases shared the same underlying TDB
                # phase (e.g. phase_map padding duplicating FCC_A1).
                # That produced an incomplete NUMPHASES x NUMPHASES table,
                # and the compiled solver segfaults reading missing
                # entries. Always emit a full square table instead.

                pair_data = eq_compositions.get((tp_i, tp_j), {})
                comp_vals = pair_data.get(tp_i, fallback_comp)
                comp_vals = [
                    (v if v is not None else fallback_comp[k])
                    for k, v in enumerate(comp_vals)
                ]

                ceq_entries.append([i, j] + comp_vals)
                cfill_entries.append([i, j] + list(comp_vals))
                cguess_entries.append(
                    [i, j] + [round(cv * 0.95, 10) for cv in comp_vals]
                )

        if failed_pairs:
            raise RuntimeError(
                f"Equilibrium calculation failed for {len(failed_pairs)} phase pair(s): "
                f"{failed_pairs}. The TDB does not define real thermodynamic data for "
                f"these phases — refusing to substitute fabricated placeholder compositions."
            )
        print(f"  [EquilibriumSolver] Generated {len(ceq_entries)} ceq, "
              f"{len(cfill_entries)} cfill, {len(cguess_entries)} c_guess "
              f"entries")
        return {
            'ceq':     ceq_entries,
            'cfill':   cfill_entries,
            'c_guess': cguess_entries,
        }
# ============================================================================
# INTERACTIVE THERMODYNAMIC DATA COLLECTION
# ============================================================================

def find_missing_phase_coverage(config, num_phases):
    """Return (covered_count, missing_pairs) -- how many phases the preset's
    existing ceq/cfill/c_guess data already covers, and which (i, j) index
    pairs (out of the full num_phases x num_phases table) still need data."""
    covered_count = 0
    for key in ('ceq', 'cfill', 'c_guess'):
        entries = config.get(key)
        if entries and isinstance(entries, list) and isinstance(entries[0], list):
            covered = set()
            for e in entries:
                covered.add(int(e[0]))
                covered.add(int(e[1]))
            covered_count = max(covered_count, max(covered) + 1)
    if covered_count >= num_phases:
        return covered_count, []
    missing_pairs = [
        (i, j)
        for i in range(num_phases)
        for j in range(num_phases)
        if i >= covered_count or j >= covered_count
    ]
    return covered_count, missing_pairs


def prompt_composition(label, n_values):
    """Ask for n_values mole fractions (0-1), reprompting until valid."""
    while True:
        raw = input(f"    {label} ({n_values} value(s), comma-separated, 0-1): ").strip()
        parts = [p.strip() for p in raw.split(",") if p.strip()]
        if len(parts) != n_values:
            print(f"      Need exactly {n_values} value(s), got {len(parts)}. Try again.")
            continue
        try:
            vals = [float(p) for p in parts]
        except ValueError:
            print("      Not all values are numbers. Try again.")
            continue
        if not all(0.0 <= v <= 1.0 for v in vals):
            print("      All values must be mole fractions between 0 and 1. Try again.")
            continue
        return vals


def prompt_for_missing_thermo_data(config, num_phases, missing_pairs):
    """Interactively ask the user for ceq/cfill/c_guess rows for every
    (i, j) pair the preset doesn't already cover. Never fabricates values."""
    n_comp = int(config.get('NUMCOMPONENTS', 2))
    n_values = max(1, n_comp - 1)
    phase_names = config.get('PHASES', [])
    if not isinstance(phase_names, list):
        phase_names = [phase_names]

    def _pname(idx):
        return phase_names[idx] if idx < len(phase_names) else f"phase {idx}"

    print()
    print("  " + "=" * 56)
    print(f"  Thermodynamic data needed for {len(missing_pairs)} phase pair(s)")
    print("  " + "=" * 56)
    print(f"  Enter equilibrium mole fraction(s) for a {n_comp}-component "
          f"system ({n_values} independent value(s) per entry).")

    new_ceq, new_cfill, new_cguess = [], [], []
    for (i, j) in missing_pairs:
        print(f"\n  Pair ({i}, {j})  [{_pname(i)} / {_pname(j)}]")
        ceq_vals = prompt_composition("ceq", n_values)

        same = input("    Use these same values for cfill and c_guess "
                      "(c_guess = ceq x 0.95)? [Y/n] ").strip().lower()
        if same in ("", "y", "yes"):
            cfill_vals = list(ceq_vals)
            cguess_vals = [round(v * 0.95, 10) for v in ceq_vals]
        else:
            cfill_vals = prompt_composition("cfill", n_values)
            cguess_vals = prompt_composition("c_guess", n_values)

        new_ceq.append([i, j] + ceq_vals)
        new_cfill.append([i, j] + cfill_vals)
        new_cguess.append([i, j] + cguess_vals)

    return {'ceq': new_ceq, 'cfill': new_cfill, 'c_guess': new_cguess}
# ============================================================================
# CONFIGURATION SCHEMA
# ============================================================================

MULTI_LINE_KEYS = frozenset({
    "DIFFUSIVITY", "EIGEN_STRAIN", "VOIGT_ISOTROPIC",
    "BOUNDARY", "BOUNDARY_VALUE",
    "ceq", "cfill", "c_guess", "slopes", "A",
    "FILLCYLINDER", "FILLSPHERE", "FILLCUBE", "FILLRECT",
    "FILLCUBERANDOM", "FILLCYLINDERRANDOM", "FILLELLIPSE",
    "FILLVORONOI2D", "FILLVORONOI3D", "FILLCUBEPATTERN",
})
BOUNDARY_FIELDS = ("phi", "mu", "c", "T")
SECTIONS = [
    ("Geometrical dimensions of the simulation domain",
     ["DIMENSION", "MESH_X", "MESH_Y", "MESH_Z"]),
    ("Discretization, space and time",
     ["DELTA_X", "DELTA_Y", "DELTA_Z", "DELTA_t"]),
    ("Number of phases and components",
     ["NUMPHASES", "NUMCOMPONENTS"]),
    ("Running and saving information",
     ["NTIMESTEPS", "NSMOOTH", "SAVET", "RESTART", "STARTTIME",
      "numworkers", "TRACK_PROGRESS"]),
    ("Component and Phase names",
     ["COMPONENTS", "PHASES"]),
    ("Material properties",
     ["GAMMA", "R", "V", "ELASTICITY", "DIFFUSIVITY",
      "EIGEN_STRAIN", "VOIGT_ISOTROPIC",
      "rho", "damping_factor", "max_iterations"]),
    ("Boundary conditions",
     ["BOUNDARY", "BOUNDARY_VALUE"]),
    ("Model-specific parameters",
     ["ISOTHERMAL", "BINARY", "TERNARY", "DILUTE", "GRAIN_GROWTH",
      "T", "WRITEFORMAT", "WRITEHDF5",
      "epsilon", "tau", "Tau"]),
    ("Anisotropy",
     ["Function_anisotropy", "Anisotropy_type", "dab", "fab",
      "Rotation_matrix"]),
    ("Potential function",
     ["Function_W", "Gamma_abc"]),
    ("Domain shifting",
     ["Shift", "Shiftj"]),
    ("Composition output",
     ["Writecomposition"]),
    ("Noise",
     ["Noise_phasefield", "Amp_Noise_Phase"]),
    ("Temperature",
     ["Equilibrium_temperature", "Filling_temperature", "Tempgrady"]),
    ("Free energy function",
     ["Function_F", "A", "ceq", "cfill", "c_guess", "slopes"]),
    ("Thermodynamic database",
     ["num_thermo_phases", "tdbfname", "tdb_phases", "phase_map"]),
    ("Elasticity experiment",
     ["eidt_mode", "comp_ff", "rad_coeff", "step_coeff",
      "step_zero", "comp_ff_rate"]),
    ("LBM",
     ["LBM", "LBM_SAVE_FREQ", "LBM_RESTART", "NU_LBM", "nu",
      "rho_LBM", "beta_c", "gy", "W0", "dt"]),
]

PARAM_ALIASES = {
    "equilibrium temperature": "Equilibrium_temperature",
    "filling temperature": "Filling_temperature",
    "number of components": "NUMCOMPONENTS",
    "number of phases": "NUMPHASES",
    "interface width": "epsilon",
    "noise amplitude": "Amp_Noise_Phase",
    "save interval": "SAVET",
    "num components": "NUMCOMPONENTS",
    "num workers": "numworkers",
    "molar volume": "V",
    "gas constant": "R",
    "num phases": "NUMPHASES",
    "time steps": "NTIMESTEPS",
    "eq temp": "Equilibrium_temperature",
    "fill temp": "Filling_temperature",
    "delta t": "DELTA_t",
    "delta x": "DELTA_X",
    "delta y": "DELTA_Y",
    "delta z": "DELTA_Z",
    "mesh x": "MESH_X",
    "mesh y": "MESH_Y",
    "mesh z": "MESH_Z",
    "tdb file": "tdbfname",
    "temperature": "T",
    "components": "NUMCOMPONENTS",
    "timesteps": "NTIMESTEPS",
    "dimension": "DIMENSION",
    "workers": "numworkers",
    "epsilon": "epsilon",
    "phases": "NUMPHASES",
    "temp": "T",
    "dim": "DIMENSION",
    "dt": "DELTA_t",
    "dx": "DELTA_X",
    "dy": "DELTA_Y",
    "dz": "DELTA_Z",
    "gamma": "GAMMA",
    "mesh": "MESH_X",
    "savet": "SAVET",
    "restart": "RESTART",
    "smooth": "NSMOOTH",
    "anisotropy": "Function_anisotropy",
    "shift": "Shift",
    "noise": "Noise_phasefield",
    "damping": "damping_factor",
    "elasticity": "ELASTICITY",
    "track progress": "TRACK_PROGRESS",
    "track_progress": "TRACK_PROGRESS",
}

# ============================================================================
# TEMPLATE ENGINE -- generate / parse .in files
# ============================================================================

def _fmt(value):
    if isinstance(value, float):
        return f"{value:g}"
    return str(value)


def _fmt_list(values):
    return "{" + ", ".join(_fmt(v) for v in values) + "}"


def _is_repeated(key, value):
    return (isinstance(value, list) and value
            and isinstance(value[0], (list, tuple)))


def generate_input_content(config):
    lines = []
    written = set()

    for header, keys in SECTIONS:
        section_lines = []
        for key in keys:
            if key not in config or key in written:
                continue
            written.add(key)
            value = config[key]
            if _is_repeated(key, value):
                for entry in value:
                    section_lines.append(f"{key} = {_fmt_list(entry)};")
            elif isinstance(value, list):
                section_lines.append(f"{key} = {_fmt_list(value)};")
            else:
                section_lines.append(f"{key} = {_fmt(value)};")
        if section_lines:
            lines.append(f"##{header}")
            lines.extend(section_lines)

    extra = []
    for key, value in config.items():
        if key in written or key.startswith("_"):
            continue
        if _is_repeated(key, value):
            for entry in value:
                extra.append(f"{key} = {_fmt_list(entry)};")
        elif isinstance(value, list):
            extra.append(f"{key} = {_fmt_list(value)};")
        else:
            extra.append(f"{key} = {_fmt(value)};")
    if extra:
        lines.append("##Additional parameters")
        lines.extend(extra)

    return "\n".join(lines) + "\n"


def generate_filling_content(filling_config):
    if not filling_config:
        return None
    lines = []
    for key, value in filling_config.items():
        if key.startswith("_"):
            continue
        if _is_repeated(key, value):
            for entry in value:
                lines.append(f"{key} = {_fmt_list(entry)};")
        elif isinstance(value, list):
            lines.append(f"{key} = {_fmt_list(value)};")
        else:
            lines.append(f"{key} = {_fmt(value)};")
    return "\n".join(lines) + "\n"

def generate_default_filling(config):
    """Build a sensible default filling config when a preset ships without
    its own Filling*.in file. Seeds phase 0 as a centered cylinder scaled
    to the mesh size."""
    mesh_x = int(config.get("MESH_X", 100))
    mesh_y = int(config.get("MESH_Y", 100))
    radius = max(3, min(mesh_x, mesh_y) // 10)
    return {
        "FILLCYLINDER": [0, mesh_x // 2, mesh_y // 2, 0, 0, radius]
    }
def rescale_filling_geometry(filling_config, base_config, config):
    """Rescale FILL* center coordinates when mesh size is overridden,
    so seed geometry stays proportionally centered in the new domain
    instead of landing outside a shrunken mesh."""
    if not filling_config:
        return filling_config

    orig_x = int(base_config.get("MESH_X", 100))
    orig_y = int(base_config.get("MESH_Y", 100))
    new_x = int(config.get("MESH_X", orig_x))
    new_y = int(config.get("MESH_Y", orig_y))

    if orig_x == new_x and orig_y == new_y:
        return filling_config  # no mesh override, nothing to rescale

    scale_x = new_x / orig_x if orig_x else 1
    scale_y = new_y / orig_y if orig_y else 1

    filling_config = copy.deepcopy(filling_config)

    # These keys have an explicit x/y center at entry[1]/entry[2] -> safe to rescale
    CENTER_BASED_KEYS = ("FILLSPHERE", "FILLCYLINDER", "FILLCYLINDERNEXLP")

    # These place particles randomly and have NO center:
    #   FILLCYLINDERRANDOM = {phase, r, vf, shield, spread}
    # entry[1] is radius, entry[2] is volume fraction -- not a y-center.
    # Treating them as coordinates corrupts the radius and can round vf
    # down to 0, killing particle placement entirely. vf is already
    # resolution-independent, so leave these untouched.
    RANDOM_FILL_KEYS = ("FILLCYLINDERRANDOM", "FILLCYLINDERRANDOMNEXLP")

    scale_r = (scale_x + scale_y) / 2  # radius scales with the average of x/y ratios

    for key in CENTER_BASED_KEYS:
        if key in filling_config:
            entries = filling_config[key]
            if entries and isinstance(entries[0], list):
                entry_list = entries
            else:
                entry_list = [entries]

            for entry in entry_list:
                if isinstance(entry, list) and len(entry) >= 3:
                    entry[1] = round(entry[1] * scale_x)  # x_center
                    entry[2] = round(entry[2] * scale_y)  # y_center
                    msg = f"  [Auto-Correction] {key} center rescaled for mesh change ({orig_x}x{orig_y} -> {new_x}x{new_y})"
                    if len(entry) >= 6:
                        old_r = entry[5]
                        entry[5] = max(1, round(entry[5] * scale_r))  # radius, floored at 1
                        msg += f", radius {old_r} -> {entry[5]}"
                    print(msg)

    for key in RANDOM_FILL_KEYS:
        if key in filling_config:
            area_ratio = scale_x * scale_y
            if area_ratio < 0.5:
                print(f"  [Note] {key} uses a resolution-independent volume fraction, so its "
                      f"parameters are left unchanged for the mesh change "
                      f"({orig_x}x{orig_y} -> {new_x}x{new_y}). If the domain shrank a lot, "
                      f"double-check that the configured volume fraction still yields at "
                      f"least 1 particle at the new size.")

    return filling_config

def parse_input_file(filepath):
    config = {}
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = re.match(r"(\S+)\s*=\s*(.+?)\s*;", line)
            if not m:
                continue
            key, raw = m.group(1), m.group(2).strip()
            value = _parse_raw(raw)
            if key in MULTI_LINE_KEYS:
                config.setdefault(key, [])
                config[key].append(value if isinstance(value, list) else [value])
            else:
                config[key] = value
    return config


def _parse_raw(raw):
    if raw.startswith("{") and raw.endswith("}"):
        inner = raw[1:-1].strip()
        if not inner:
            return []
        parts = re.split(r"[,\t]+", inner)
        return [_coerce(p.strip()) for p in parts if p.strip()]
    return _coerce(raw)


def _coerce(s):
    try:
        if "." not in s and "e" not in s.lower():
            return int(s)
        return float(s)
    except (ValueError, AttributeError):
        return s

# ============================================================================
# AUTO-DISCOVERY -- scan Example_Systems for alloy presets
# ============================================================================

def discover_systems(base_dir):
    """Scan Example_Systems/ and build configs from the actual .in files on disk."""
    systems = {}
    example_dir = os.path.join(base_dir, "Example_Systems")
    if not os.path.isdir(example_dir):
        return systems

    for name in sorted(os.listdir(example_dir)):
        sys_dir = os.path.join(example_dir, name)
        if not os.path.isdir(sys_dir):
            continue

        input_files = sorted(glob.glob(os.path.join(sys_dir, "Input*.in")))
        if not input_files:
            continue

        primary_input = input_files[0]
        try:
            config = parse_input_file(primary_input)
        except Exception:
            continue

        if not config:
            continue

        filling_files = sorted(glob.glob(os.path.join(sys_dir, "Fill*.in")))
        filling_config = None
        if filling_files:
            try:
                filling_config = parse_input_file(filling_files[0])
            except Exception:
                pass

        if not filling_config:
            filling_config = generate_default_filling(config)

        components = config.get("COMPONENTS", [])
        phases = config.get("PHASES", [])

        systems[name] = {
            "config": config,
            "filling": filling_config,
            "input_file": os.path.basename(primary_input),
            "filling_file": os.path.basename(filling_files[0]) if filling_files else None,
            "all_inputs": [os.path.basename(f) for f in input_files],
            "all_fillings": [os.path.basename(f) for f in filling_files],
            "dir": sys_dir,
            "components": components if isinstance(components, list) else [components],
            "phases": phases if isinstance(phases, list) else [phases],
        }

    return systems

# ============================================================================
# DISPLAY HELPERS
# ============================================================================

def show_systems(systems):
    print()
    print("=" * 60)
    print("  AVAILABLE ALLOY SYSTEMS")
    print("=" * 60)
    for i, (name, info) in enumerate(systems.items(), 1):
        comps = ", ".join(str(c) for c in info["components"])
        n_phases = info["config"].get("NUMPHASES", "?")
        temp = info["config"].get("T", "?")
        eq_temp = info["config"].get("Equilibrium_temperature", "?")
        print(f"  {i:2d}. {name}")
        print(f"      Components: [{comps}]  Phases: {n_phases}  T: {temp}  Eq.T: {eq_temp}")
    print()


def show_config_diff(base_config, overrides):
    if not overrides:
        print("  No parameter changes.")
        return
    print("  Parameter changes:")
    for k, v in sorted(overrides.items()):
        old = base_config.get(k, "(unset)")
        if isinstance(old, list):
            old = _fmt_list(old)
        print(f"    {k}: {old} -> {v}")


def show_key_params(config):
    print("  Key parameters of the final config:")
    show_keys = [
        "DIMENSION", "MESH_X", "MESH_Y", "MESH_Z",
        "DELTA_X", "DELTA_Y", "DELTA_Z", "DELTA_t",
        "NUMPHASES", "NUMCOMPONENTS",
        "NTIMESTEPS", "SAVET",
        "T", "Equilibrium_temperature", "Filling_temperature",
        "epsilon", "Function_F",
    ]
    for k in show_keys:
        if k in config:
            v = config[k]
            if isinstance(v, list):
                v = _fmt_list(v)
            print(f"    {k} = {v}")
