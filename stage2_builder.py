#!/usr/bin/env python3
"""
stage2_builder.py
===================
Takes the structured spec produced by stage1_translator.py and turns it
into real Input.in / Filling.in files microsim_gp can run.

This stage does NOT invent thermodynamic data. It always starts from the
existing preset values under Example_Systems/ (ceq/cfill/c_guess, TDB
choices, etc.) and only reaches for EquilibriumSolver (pycalphad) when
the user asked for more phases than the preset's existing data covers --
i.e. only when *more* is genuinely required, exactly as requested.

Usage:
  python3 stage2_builder.py spec.json -o build_result.json

Env vars:
  MICROSIM_DIR   Path to microsim_gp directory (default: /app/microsim_gp)
"""

import os
import sys
import json
import math
import copy
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from shared_core import (
    EquilibriumSolver, MULTI_LINE_KEYS, discover_systems,
    generate_input_content, generate_filling_content,
    rescale_filling_geometry, show_config_diff, show_key_params,
    find_missing_phase_coverage, prompt_for_missing_thermo_data,
)

def build_config(base_config, overrides=None, interactive=True, allow_solver=False):
    """Layer user overrides and auto-pad all known phase-dependent arrays."""
    boundary_field_overrides = None
    if overrides:
        boundary_field_overrides = overrides.pop("_boundary_field_overrides", None)
        for bkey in ("BOUNDARY", "BOUNDARY_VALUE"):
            if overrides.get(bkey) == "__STRUCTURED_PENDING__":
                overrides.pop(bkey)

    config = copy.deepcopy(base_config)
    if overrides:
        for key, val in list(overrides.items()):
            existing = base_config.get(key)
            # A MULTI_LINE_KEYS parameter (DIFFUSIVITY, EIGEN_STRAIN, ceq,
            # etc.) is stored as a list of structured per-phase entries,
            # e.g. [diagonal_flag, phase_index, D11, D22, ...]. If the
            # preset already stores it that way but an override tries to
            # replace it with a bare scalar or flat list, silently
            # accepting that would corrupt the structure the compiled
            # solver expects -- exactly the kind of malformed input that
            # crashes microsim_gp instead of failing gracefully. Reject
            # it here instead of letting it reach the padding logic.
            if (key in MULTI_LINE_KEYS and isinstance(existing, list)
                    and existing and isinstance(existing[0], list)
                    and not (isinstance(val, list) and val and isinstance(val[0], list))):
                print(f"  [Auto-Correction] Ignoring override for '{key}' -- this preset "
                      f"stores it as structured per-phase data (e.g. "
                      f"{existing[0]}), but the override was a plain value "
                      f"({val!r}). Refusing to flatten it, since that would "
                      f"produce a config the solver can't safely parse. "
                      f"Keeping the preset's existing '{key}' data unchanged.")
                overrides.pop(key)
        config.update(overrides)

    if boundary_field_overrides:
        for bkey, field_map in boundary_field_overrides.items():
            existing_rows = config.get(bkey)
            if not isinstance(existing_rows, list):
                existing_rows = []
            rows_by_field = {row[0]: row for row in existing_rows if isinstance(row, list) and row}
            for field, vals in field_map.items():
                rows_by_field[field] = [field] + vals
            ordered = []
            for f in BOUNDARY_FIELDS:
                if f in rows_by_field:
                    ordered.append(rows_by_field.pop(f))
            ordered.extend(rows_by_field.values())  # preserve any unexpected extra rows
            config[bkey] = ordered
            print(f"  [Auto-Correction] {bkey} updated for field(s): {', '.join(field_map.keys())} (other fields kept from preset)")

    # 1. Configuration Registry
    num_phases = int(config.get('NUMPHASES', 1))
    
    # Comprehensive list of parameters that usually require phase-specific values
    # Added 'phase_map' to the list below
    num_pairs = max(1, num_phases * (num_phases - 1) // 2)

    # Per-phase parameters (one value per phase). NOTE: epsilon and tau are
    # deliberately excluded -- they are always scalar values in microsim_gp,
    # never per-phase arrays. Including them here causes the compiled
    # solver's own validator to reject the generated Input.in with
    # "epsilon must be a number, got: list".
    PARAMS_TO_PAD = [
        'ELASTICITY', 'DIFFUSIVITY', 'EIGEN_STRAIN', 'VOIGT_ISOTROPIC',
        'rho', 'damping_factor', 'A', 'ceq', 'cfill', 'c_guess', 'slopes', 'phase_map'
    ]

    # Per-phase-pair parameters (e.g. surface energies)
    PAIRWISE_PARAMS = ['GAMMA', 'Tau', 'dab', 'fab']

    PHASE_INDEX_POSITION = {
        'DIFFUSIVITY': 1,
    }
    # Gamma_abc is a three-phase-junction term: one value per unique triplet
    # of phases (analogous to GAMMA being one value per pair). Required count
    # is C(NUMPHASES, 3), floored at 1 -- confirmed against the binary's own
    # error messages: NUMPHASES=3 -> 1 required, NUMPHASES=4 -> 4 required.
    if 'Gamma_abc' in config:
        val = config['Gamma_abc']
        if not isinstance(val, list):
            val = [val]
        num_triplets = max(1, math.comb(num_phases, 3))
        if len(val) < num_triplets:
            fill_val = val[0] if val else 0.1
            val = val + [fill_val] * (num_triplets - len(val))
            print(f"  [Auto-Correction] Gamma_abc too short. Padding to {num_triplets} (per phase-triplet) with {fill_val}")
        elif len(val) > num_triplets:
            print(f"  [Auto-Correction] Gamma_abc too long ({len(val)}). Trimming to {num_triplets} (per phase-triplet)")
            val = val[:num_triplets]
        config['Gamma_abc'] = val

    # 2. Dynamic Auto-Padding — per-phase


    # 2. Dynamic Auto-Padding — per-phase
    for key in PARAMS_TO_PAD:
        if key in config:
            val = config[key]
            if not isinstance(val, list):
                val = [val]
            is_pair_indexed = key in ('ceq', 'cfill', 'c_guess')
            if len(val) < num_phases and not is_pair_indexed:
                if key in MULTI_LINE_KEYS and val and isinstance(val[0], list):
                    idx_pos = PHASE_INDEX_POSITION.get(key, 0)
                    template = val[0]
                    print(f"  [Auto-Correction] {key} too short ({len(val)}). Padding to {num_phases} (indexed, cloned from phase {template[idx_pos]}, inserted before final/reference phase)")
                    while len(val) < num_phases:
                        insert_pos = len(val) - 1
                        new_entry = list(template)
                        new_entry[idx_pos] = insert_pos
                        val.insert(insert_pos, new_entry)
                else:
                    fill_val = val[0] if len(val) > 0 else "None"
                    print(f"  [Auto-Correction] {key} too short ({len(val)}). Padding to {num_phases} with {fill_val}, inserted before final/reference phase")
                    while len(val) < num_phases:
                        val.insert(len(val) - 1, fill_val)
                config[key] = val
            elif val and isinstance(val[0], list):
                if is_pair_indexed:
                    # ceq/cfill/c_guess entries reference two phase indices [a, b, ...];
                    # drop any entry that refers to a phase index no longer in range.
                    kept = [e for e in val if int(e[0]) < num_phases and int(e[1]) < num_phases]
                    if len(kept) != len(val):
                        print(f"  [Auto-Correction] {key} trimmed: removed {len(val)-len(kept)} entries referencing phases >= {num_phases}")
                        config[key] = kept
                elif len(val) > num_phases:
                    # Indexed single-phase entries; keep only those whose
                    # phase index (position varies by key) is still valid.
                    idx_pos = PHASE_INDEX_POSITION.get(key, 0)
                    kept = [e for e in val if int(e[idx_pos]) < num_phases]
                    print(f"  [Auto-Correction] {key} too long ({len(val)}). Trimming to {len(kept)} (indexed)")
                    config[key] = kept
            elif len(val) > num_phases:
                print(f"  [Auto-Correction] {key} too long ({len(val)}). Trimming to {num_phases}")
                config[key] = val[:num_phases]

    # 2b. Dynamic Auto-Padding — pairwise
    # Safety net: epsilon/tau must always be scalars in microsim_gp. If
    # they arrived as a list from anywhere (an old spec.json, a stray
    # override), collapse back to a single value instead of writing an
    # invalid Input.in.
    for key in ("epsilon", "tau"):
        if key in config and isinstance(config[key], list):
            if config[key]:
                config[key] = config[key][0]
                print(f"  [Auto-Correction] {key} must be scalar; using first value.")
            else:
                del config[key]

    for key in PAIRWISE_PARAMS:
        if key in config:
            val = config[key]
            if not isinstance(val, list):
                val = [val]
            if len(val) < num_pairs:
                fill_val = val[0] if len(val) > 0 else 0.1
                print(f"  [Auto-Correction] {key} too short ({len(val)}). Padding to {num_pairs} (pairwise) with {fill_val}")
                while len(val) < num_pairs:
                    val.append(fill_val)
                config[key] = val
            elif len(val) > num_pairs:
                print(f"  [Auto-Correction] {key} too long ({len(val)}). Trimming to {num_pairs} (pairwise)")
                config[key] = val[:num_pairs]

    # 3. Handle PHASES naming specifically
    if 'PHASES' in config:
        phases = config['PHASES']
        if not isinstance(phases, list):
            phases = [phases]
        if len(phases) < num_phases:
            missing = num_phases - len(phases)
            new_names = [f"beta_{i+1}" for i in range(missing)]
            print(f"  [Auto-Correction] Adding {missing} default phase(s): {', '.join(new_names)} (inserted before the final/reference phase)")
            # Insert before the last phase so whatever was originally last
            # (assumed liquid/reference by the compiled solver) stays last.
            insert_at = len(phases) - 1
            phases[insert_at:insert_at] = new_names
            config['PHASES'] = phases
        elif len(phases) > num_phases:
            print(f"  [Auto-Correction] PHASES too long ({len(phases)}). Trimming to {num_phases}")
            config['PHASES'] = phases[:num_phases]
    
    # 4. Dynamic solver — extend ceq/cfill/c_guess when phases exceed
    #    the pre-defined template data.  The original hard error is kept
    #    as a fallback if the solver is unavailable or fails.
    covered_count, missing_pairs = find_missing_phase_coverage(config, num_phases)

    if missing_pairs:
        if interactive:
            print(f"\n  [builder] This preset only has thermodynamic data for "
                  f"{covered_count} phase(s); NUMPHASES={num_phases} needs more.")
            try:
                computed = prompt_for_missing_thermo_data(config, num_phases, missing_pairs)
            except (EOFError, KeyboardInterrupt):
                raise ValueError(
                    f"Thermodynamic data entry cancelled -- NUMPHASES={num_phases} "
                    f"needs ceq/cfill/c_guess data for phase(s) "
                    f"{covered_count}..{num_phases - 1} that this preset doesn't have."
                )
            for key in ('ceq', 'cfill', 'c_guess'):
                config.setdefault(key, [])
                config[key] = config[key] + computed[key]
        elif allow_solver:
            try:
                _solver = EquilibriumSolver(
                    base_dir=os.environ.get("MICROSIM_DIR", "/app/microsim_gp")
                )
                _computed = _solver.compute_params(config, num_phases)
                for _ckey in ('ceq', 'cfill', 'c_guess'):
                    if _ckey in _computed:
                        config[_ckey] = _computed[_ckey]
            except ImportError as e:
                print(f"  [EquilibriumSolver] Unavailable ({e}).")
            except Exception as e:
                print(f"  [EquilibriumSolver] Computation failed: {e}.")
        # else: leave it missing -- the validation loop below raises a
        # clear error naming exactly which phase(s) lack data.

    # Original validation — catches cases where solver was skipped or failed
    for key in ('ceq', 'cfill', 'c_guess'):
        if key in config:
            entries = config[key]
            if entries and isinstance(entries[0], list):
                covered = set()
                for entry in entries:
                    covered.add(int(entry[0]))
                    covered.add(int(entry[1]))
                max_covered = max(covered) + 1
                if num_phases > max_covered:
                    print(f"  [ERROR] This preset's '{key}' data only covers {max_covered} phase(s).")
                    print(f"  [ERROR] Increasing NUMPHASES to {num_phases} would leave phase(s) "
                          f"{max_covered}..{num_phases-1} without real composition data, "
                          f"which would likely crash the solver. Aborting before running.")
                    raise ValueError(f"Insufficient '{key}' data for NUMPHASES={num_phases}")

    return config

def write_simulation_files(config, filling_config=None,
                           output_dir="/app/microsim_gp",
                           input_name="Input_generated.in",
                           filling_name="Filling_generated.in"):
    os.makedirs(output_dir, exist_ok=True)

    input_path = os.path.join(output_dir, input_name)
    with open(input_path, "w") as f:
        f.write(generate_input_content(config))

    filling_path = None
    if filling_config:
        filling_path = os.path.join(output_dir, filling_name)
        with open(filling_path, "w") as f:
            f.write(generate_filling_content(filling_config))

    return input_path, filling_path


# ============================================================================
# MAIN ENTRY POINT
# ============================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("spec", help="Path to spec.json produced by stage1_translator.py")
    ap.add_argument("-o", "--output", default="build_result.json",
                     help="Where to write the build result JSON")
    ap.add_argument("--base-dir", default=os.environ.get("MICROSIM_DIR", "/app/microsim_gp"))
    ap.add_argument("--non-interactive", action="store_true",
                     help="Never prompt for missing thermodynamic data; use --auto-solve or fail")
    ap.add_argument("--auto-solve", action="store_true",
                     help="If data is missing and prompting is off, fall back to "
                          "automatic pycalphad computation instead of failing")
    args = ap.parse_args()

    with open(args.spec) as f:
        spec = json.load(f)

    system_name = spec.get("system")
    overrides = dict(spec.get("overrides") or {})

    systems = discover_systems(args.base_dir)
    if not systems:
        print(f"  [builder] ERROR: no systems found under "
              f"{os.path.join(args.base_dir, 'Example_Systems')}")
        sys.exit(1)

    if system_name and system_name in systems:
        base_config = systems[system_name]["config"]
        filling = systems[system_name].get("filling")
    else:
        if system_name:
            print(f"  [builder] System '{system_name}' not recognized. "
                  f"Available: {', '.join(systems.keys())}")
            sys.exit(1)
        system_name = next(iter(systems))
        print(f"  [builder] No system specified in spec -- defaulting to '{system_name}'")
        base_config = systems[system_name]["config"]
        filling = systems[system_name].get("filling")

    print(f"  [builder] System: {system_name}")
    show_config_diff(base_config, overrides)

    try:
        config = build_config(base_config, overrides,
                               interactive=not args.non_interactive,
                               allow_solver=args.auto_solve)
    except ValueError as e:
        print(f"\n  [builder] ERROR: cannot build this configuration: {e}")
        sys.exit(1)

    filling = rescale_filling_geometry(filling, base_config, config)

    show_key_params(config)

    input_path, filling_path = write_simulation_files(
        config, filling, output_dir=args.base_dir,
    )
    print(f"  [builder] Wrote: {input_path}")
    if filling_path:
        print(f"  [builder] Wrote: {filling_path}")

    result = {
        "system": system_name,
        "base_dir": args.base_dir,
        "input_path": input_path,
        "filling_path": filling_path,
        "config": config,
        "filling": filling,
        "dimension": config.get("DIMENSION", 2),
    }
    with open(args.output, "w") as f:
        json.dump(result, f, indent=2)
    print(f"  [builder] Wrote build result to {args.output}")


if __name__ == "__main__":
    main()
