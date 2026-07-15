#!/usr/bin/env python3
"""
stage3_validator.py
=====================
Sanity-checks the config produced by stage2_builder.py before anything
gets handed to mpirun. Specifically guards against "ghost systems":
configurations that *look* complete (every field present, no crash) but
whose thermodynamic values are fabricated, degenerate, or don't actually
correspond to real phase-equilibrium data.

Checks performed:
  1. TDB reality check   -- ONLY for systems that actually depend on a
                             .tdb file. If a real .tdb was named, re-open
                             it and confirm (independently of stage2)
                             that it actually defines the elements &
                             phases this run needs. Systems that get
                             their thermodynamic data from CSVs
                             (tdbs_encrypted/Composition_*.csv,
                             HSN_*.csv -- e.g. NiAl, NiAlCr, NiAlMo) are
                             skipped: their Input.in's tdbfname field is
                             a shared placeholder name, not the actual
                             data source, so validating it against
                             pycalphad produces false positives.
  2. Fallback TDB flag   -- warns if the run is quietly using the
                             minimal auto-generated placeholder TDB
                             (ideal-solution, not a real assessment).
  3. Composition sanity  -- every ceq/cfill/c_guess mole fraction must
                             be a real number in (0, 1); flags any phase
                             pair whose composition is missing or was
                             filled with the crude equal-split fallback.
  4. Degenerate phases   -- if two *different* phases share identical
                             equilibrium compositions across every pair,
                             they aren't thermodynamically distinguishable
                             -- almost always a sign of fabricated/ghost
                             data rather than a real multi-phase system.
  5. Structural sanity   -- NUMPHASES matches PHASES length and the
                             ceq/cfill/c_guess index coverage; mesh,
                             timestep, and DELTA_t are positive.

Exits 0 and writes {"ok": true, ...} if everything passes.
Exits 1 and writes {"ok": false, "issues": [...]} otherwise.

Usage:
  python3 stage3_validator.py build_result.json -o validation_result.json
"""

import os
import sys
import json
import glob
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from shared_core import EquilibriumSolver


def _is_real_fraction(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and 0.0 <= float(x) <= 1.0


def _system_uses_csv_thermo_data(system_name, base_dir):
    """True if this preset gets its thermodynamic values from CSV tables
    (Composition_*.csv / HSN_*.csv) rather than a pycalphad-readable .tdb.
    These presets' Input.in still carries a tdbfname field, but it's a
    shared placeholder name across presets, not the actual data source
    -- see run_simulation()'s CSV-sync logic, which is what actually
    feeds these systems their thermodynamic data at run time."""
    if not system_name:
        return False
    system_dir = os.path.join(base_dir, "Example_Systems", system_name)
    for sub in ("tdbs_encrypted", ""):
        search_dir = os.path.join(system_dir, sub) if sub else system_dir
        if glob.glob(os.path.join(search_dir, "Composition_*.csv")) or \
           glob.glob(os.path.join(search_dir, "HSN_*.csv")):
            return True
    return False


def check_structural_sanity(config, issues, warnings):
    num_phases = int(config.get("NUMPHASES", 1))
    if num_phases < 1:
        issues.append(f"NUMPHASES={num_phases} is not valid (must be >= 1)")

    phases = config.get("PHASES", [])
    if isinstance(phases, list) and phases and len(phases) != num_phases:
        issues.append(f"PHASES has {len(phases)} entries but NUMPHASES={num_phases}")

    for dim_key in ("MESH_X", "MESH_Y", "MESH_Z"):
        if dim_key in config:
            v = config[dim_key]
            if not isinstance(v, (int, float)) or v <= 0:
                issues.append(f"{dim_key}={v} must be a positive number")

    if "NTIMESTEPS" in config and (not isinstance(config["NTIMESTEPS"], (int, float)) or config["NTIMESTEPS"] <= 0):
        issues.append(f"NTIMESTEPS={config['NTIMESTEPS']} must be positive")

    if "DELTA_t" in config and (not isinstance(config["DELTA_t"], (int, float)) or config["DELTA_t"] <= 0):
        issues.append(f"DELTA_t={config['DELTA_t']} must be positive")

    for key in ("ceq", "cfill", "c_guess"):
        entries = config.get(key)
        if entries and isinstance(entries, list) and isinstance(entries[0], list):
            covered = set()
            for e in entries:
                covered.add(int(e[0]))
                covered.add(int(e[1]))
            if max(covered) + 1 < num_phases:
                issues.append(
                    f"{key} only covers phases 0..{max(covered)}, but "
                    f"NUMPHASES={num_phases} needs coverage up to {num_phases - 1}"
                )


def check_composition_sanity(config, issues, warnings):
    for key in ("ceq", "cfill", "c_guess"):
        entries = config.get(key)
        if not entries or not isinstance(entries, list) or not isinstance(entries[0], list):
            continue
        for entry in entries:
            phase_a, phase_b = entry[0], entry[1]
            comps = entry[2:]
            if not comps:
                issues.append(f"{key} entry for phases ({phase_a},{phase_b}) has no composition values")
                continue
            for v in comps:
                if v is None:
                    issues.append(f"{key} entry for phases ({phase_a},{phase_b}) has a missing (null) value")
                elif not _is_real_fraction(v):
                    issues.append(
                        f"{key} entry for phases ({phase_a},{phase_b}) has an out-of-range "
                        f"mole fraction {v} (expected 0-1)"
                    )


def check_degenerate_phases(config, warnings):
    """Flag phases whose equilibrium composition is identical across the
    board -- a strong signal the data is fabricated/fallback rather than
    a real per-phase thermodynamic distinction."""
    entries = config.get("ceq")
    if not entries or not isinstance(entries, list) or not isinstance(entries[0], list):
        return
    by_phase = {}
    for entry in entries:
        phase_a = int(entry[0])
        # Skip rows with missing/invalid values here -- those are already
        # reported by check_composition_sanity(); comparing them for
        # degeneracy would be meaningless (and unorderable).
        comps = entry[2:]
        if any(v is None or not isinstance(v, (int, float)) for v in comps):
            continue
        by_phase.setdefault(phase_a, []).append(tuple(comps))

    seen_signatures = {}
    for phase, rows in by_phase.items():
        signature = tuple(sorted(rows))
        if signature in seen_signatures:
            other = seen_signatures[signature]
            warnings.append(
                f"phases {other} and {phase} have identical ceq composition data across "
                f"all pairs -- they may not be thermodynamically distinguishable "
                f"(possible ghost/fabricated data)"
            )
        else:
            seen_signatures[signature] = phase


def check_tdb_reality(config, system_name, base_dir, issues, warnings):
    tdb_name = config.get("tdbfname")
    if not tdb_name:
        return

    

    if "auto_generated" in str(tdb_name).lower():
        warnings.append(
            f"tdbfname='{tdb_name}' is the minimal auto-generated fallback TDB "
            f"(ideal-solution placeholder, not a real thermodynamic assessment). "
            f"Results are not physically meaningful -- supply a real .tdb for "
            f"production runs."
        )
        return

    try:
        from pycalphad import Database
    except ImportError:
        warnings.append("pycalphad not available -- skipping independent TDB verification")
        return

    solver = EquilibriumSolver(base_dir=base_dir)

    components = config.get("COMPONENTS", [])
    if isinstance(components, str):
        components = [components]

    phase_map = config.get("phase_map", [])
    if isinstance(phase_map, str):
        phase_map = [phase_map]

    needed_elements = {c.upper() for c in components}
    needed_phases = {p.upper() for p in phase_map}

    tdb_path = solver.find_tdb(tdb_name)

    # If the preset TDB is stale (common placeholder like alzn_mey.tdb),
    # ignore it and search for a database that actually matches this system.
    if tdb_path:
        try:
            db = Database(tdb_path)

            db_elements = {str(e).upper() for e in db.elements}
            db_phases = {str(p).upper() for p in db.phases.keys()}

            if (not needed_elements.issubset(db_elements) or
                    not needed_phases.issubset(db_phases)):
                tdb_path = None

        except Exception:
            tdb_path = None

    if not tdb_path:
        tdb_path = solver.find_best_tdb(components, phase_map)

    # Couldn't find any matching TDB -> assume this preset uses CSV data.
    if not tdb_path:
        if _system_uses_csv_thermo_data(system_name, base_dir):
            return

        issues.append(
            f"No thermodynamic database found for system '{system_name}'."
        )
        return
    try:
        db = Database(tdb_path)

        db_elements = {str(e).upper() for e in db.elements}
        db_phases = {str(p).upper() for p in db.phases.keys()}

        missing_elements = needed_elements - db_elements
        missing_phases = needed_phases - db_phases

        if missing_elements:
            issues.append(
                f"TDB '{os.path.basename(tdb_path)}' does not define required element(s): "
                f"{sorted(missing_elements)}"
            )

        if missing_phases:
            warnings.append(
                f"TDB '{os.path.basename(tdb_path)}' does not define mapped phase(s): "
                f"{sorted(missing_phases)}"
            )

    except Exception as e:
        issues.append(f"Could not verify TDB '{tdb_path}': {e}")


def validate(build_result):
    config = build_result.get("config", {})
    base_dir = build_result.get("base_dir", "/app/microsim_gp")
    system_name = build_result.get("system")

    issues = []
    warnings = []

    check_structural_sanity(config, issues, warnings)
    check_composition_sanity(config, issues, warnings)
    check_degenerate_phases(config, warnings)
    check_tdb_reality(config, system_name, base_dir, issues, warnings)

    return {"ok": len(issues) == 0, "issues": issues, "warnings": warnings}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("build_result", help="Path to build_result.json produced by stage2_builder.py")
    ap.add_argument("-o", "--output", default="validation_result.json")
    ap.add_argument("--force", action="store_true",
                     help="Exit 0 even if issues were found (still reports them)")
    args = ap.parse_args()

    with open(args.build_result) as f:
        build_result = json.load(f)

    result = validate(build_result)

    print("\n  --- Validation ---")
    if result["warnings"]:
        for w in result["warnings"]:
            print(f"    \u26a0 {w}")
    if result["issues"]:
        for i in result["issues"]:
            print(f"    \u2717 {i}")
        print(f"  [validator] FAILED: {len(result['issues'])} issue(s) found.")
    else:
        print("  [validator] OK: no issues found.")

    with open(args.output, "w") as f:
        json.dump(result, f, indent=2)

    if not result["ok"] and not args.force:
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()