#!/usr/bin/env python3
"""
run_agents.py
==============
Thin orchestrator that chains the three pipeline stages together and
then runs the simulation:

  stage1_translator.py   user paragraph  -> spec.json
  stage2_builder.py       spec.json      -> Input.in / Filling.in + build_result.json
  stage3_validator.py     build_result   -> validation_result.json
  run_simulation()        (this file)    -> mpirun microsim_gp

Each stage is a standalone script -- run_agents.py just shells out to
them in order and stops the pipeline if a stage fails. This means each
stage can be edited, tested, or replaced independently.

Usage:
  Interactive:   python3 run_agents.py
  One-shot:      python3 run_agents.py "NiAl with equilibrium temperature 1600K"

Env vars:
  MICROSIM_DIR   Path to microsim_gp directory (default: /app/microsim_gp)
  MPI_NP         Number of MPI processes (default: 4)
  MPI_GRID       MPI grid decomposition (default: '2 2')
  OLLAMA_HOST    Ollama server URL (default: http://localhost:11434)
  OLLAMA_MODEL   Ollama model for translation (default: gemma2:2b)
"""

import os
import re
import sys
import json
import glob
import shutil
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ============================================================================
# SIMULATION RUNNER -- direct subprocess, no LLM code generation
# ============================================================================

def run_simulation(input_path, filling_path=None,
                   np_procs=4, grid="2 2",
                   work_dir="/app/microsim_gp",
                   system_name=None, base_dir=None,
                   dimension=2):
    """Run microsim_gp via mpirun directly. Output goes to DATA/ inside work_dir."""
    data_dir = os.path.join(work_dir, "DATA")
    os.makedirs(data_dir, exist_ok=True)

    # Sync this system's thermodynamic data tables into the path
    # the compiled solver actually reads from (relative to work_dir)
    if system_name and base_dir:
        system_dir = os.path.join(base_dir, "Example_Systems", system_name)
        system_tdbs_dir = os.path.join(system_dir, "tdbs_encrypted")
        if not os.path.isdir(system_tdbs_dir):
            # Some presets (e.g. NiAlCr) keep their CSV thermodynamic tables
            # directly in the system folder instead of a tdbs_encrypted/ subfolder.
            system_tdbs_dir = system_dir

        global_tdbs_dir = os.path.join(work_dir, "tdbs_encrypted")
        csv_files = glob.glob(os.path.join(system_tdbs_dir, "Composition_*.csv")) + \
                    glob.glob(os.path.join(system_tdbs_dir, "HSN_*.csv"))
        if csv_files:
            os.makedirs(global_tdbs_dir, exist_ok=True)
            for fpath in csv_files:
                shutil.copy2(fpath, os.path.join(global_tdbs_dir, os.path.basename(fpath)))
            print(f"  Synced tdbs_encrypted/ from {system_name} preset ({len(csv_files)} files)") 

            # --- Validate requested temperatures against the tabulated
            #     data range -----------------------------------------------
            # The compiled solver's GSL cubic splines abort (SIGABRT,
            # "interp.c: ERROR: interpolation error") if evaluated outside
            # the tabulated T range in these CSVs -- GSL refuses to
            # extrapolate. Catch this here, before invoking the binary,
            # instead of letting 4 MPI ranks crash mid-run.
            t_min, t_max = None, None
            for fpath in csv_files:
                try:
                    with open(fpath) as fh:
                        next(fh)  # skip header
                        for line in fh:
                            line = line.strip()
                            if not line:
                                continue
                            first_field = line.split(",")[0]
                            t_val = float(first_field)
                            if t_min is None or t_val < t_min:
                                t_min = t_val
                            if t_max is None or t_val > t_max:
                                t_max = t_val
                except Exception:
                    continue

            if t_min is not None and t_max is not None:
                requested_temps = {}
                for fpath_to_check in [p for p in (input_path, filling_path) if p]:
                    try:
                        with open(fpath_to_check) as fh:
                            content = fh.read()
                        for key in ("T", "Equilibrium_temperature", "Filling_temperature"):
                            m = re.search(rf"^{key}\s*=\s*([-\d.eE+]+)", content, re.MULTILINE)
                            if m:
                                requested_temps[key] = float(m.group(1))
                    except Exception:
                        continue

                out_of_range = {
                    k: v for k, v in requested_temps.items()
                    if v < t_min or v > t_max
                }
                if out_of_range:
                    print(f"\n  [ERROR] Requested temperature(s) fall outside the tabulated "
                          f"thermodynamic data range for '{system_name}' "
                          f"({t_min:.1f}K - {t_max:.1f}K):")
                    for k, v in out_of_range.items():
                        print(f"    {k} = {v} K  (outside {t_min:.1f}-{t_max:.1f} K)")
                    print(f"  The compiled solver uses spline interpolation and will abort "
                          f"if asked to evaluate outside this range -- it cannot extrapolate. "
                          f"Choose a temperature within {t_min:.1f}-{t_max:.1f} K, or supply "
                          f"thermodynamic data covering the range you need.")
                    return 1

    # --- Run the project's own built-in validators as an additional
    #     pre-flight layer, on top of the checks above -----------------
    # NOTE: validate_thermo.py is deliberately NOT run here. It has a
    # false-positive bug (requires len(tdb_phases) == len(phase_map),
    # which is wrong whenever multiple sim-phases share one thermo
    # phase -- exactly the case for every successful TDB-based run in
    # this project, e.g. NiAl's tdb_phases={FCC_A1,LIQUID} (len 2) vs
    # phase_map={FCC_A1,LIQUID,FCC_A1} (len 3)). Using it would reject
    # valid, working configurations.
    _binary_dir_for_validators = os.environ.get("MICROSIM_DIR", "/app/microsim_gp")
    _validator_dir = os.path.join(_binary_dir_for_validators, "scripts", "validator")
    if os.path.isdir(_validator_dir):
        if _validator_dir not in sys.path:
            sys.path.insert(0, _validator_dir)
        try:
            from validate_input import InputValidator
            from validate_filling import FillingValidator
        except ImportError:
            InputValidator = FillingValidator = None

        if InputValidator is not None:
            print("\n  --- Built-in validators (scripts/validator/) ---")

            input_validator = InputValidator(input_path)
            input_ok = input_validator.parse() and input_validator.validate()
            if input_ok:
                print("  [Validator] Input file: OK")
            else:
                print(f"  [Validator] Input file: {len(input_validator.errors)} error(s)")
                for err in input_validator.errors:
                    print(f"    \u2717 {err}")
            for warn in input_validator.warnings:
                print(f"    \u26a0 {warn}")

            # Filling file: advisory only. Its VALID_FILL_TYPES whitelist
            # is missing several fill types this project actually uses
            # (FILLCUBERANDOM, FILLRECT, FILLELLIPSE, FILLCUBEPATTERN),
            # so a hard block here would false-flag legitimate configs.
            if filling_path and os.path.exists(filling_path) and FillingValidator is not None:
                mesh_size = (
                    input_validator.params.get('MESH_X'),
                    input_validator.params.get('MESH_Y'),
                    input_validator.params.get('MESH_Z'),
                )
                filling_validator = FillingValidator(filling_path, mesh_size)
                if filling_validator.parse():
                    filling_validator.validate()
                    if filling_validator.errors:
                        print(f"  [Validator] Filling file: {len(filling_validator.errors)} "
                              f"issue(s) (advisory only -- fill-type whitelist may be incomplete):")
                        for err in filling_validator.errors:
                            print(f"    \u26a0 {err}")
                    else:
                        print("  [Validator] Filling file: OK")
                    for warn in filling_validator.warnings:
                        print(f"    \u26a0 {warn}")

            print("  " + "-" * 49)
            # NOTE: InputValidator is advisory-only, not hard-blocking.
            # Confirmed false positive: its REQUIRED_SCALARS table treats
            # epsilon/tau as scalars, but this project always generates
            # them as per-phase arrays (epsilon = {v1, v2, v3};) for any
            # NUMPHASES > 1 -- which would hard-block nearly every real
            # run in this project, not just edge cases.
            if not input_ok:
                print("\n  [Validator] Input file has advisory findings above "
                      "(some checks in this validator assume a different "
                      "config convention than this project uses -- e.g. "
                      "epsilon/tau as scalars vs. per-phase arrays here). "
                      "Not blocking the run.")

            

    binary_dir = os.environ.get("MICROSIM_DIR", "/app/microsim_gp")
    binary = os.path.join(binary_dir, "microsim_gp")
    if not os.path.isfile(binary):
        print(f"ERROR: Simulation binary not found at {binary}")
        print(f"Contents of {binary_dir}:")
        try:
            for f in os.listdir(binary_dir):
                print(f"  {f}")
        except Exception as e:
            print(f"  (could not list directory: {e})")
        return 1

    input_basename = os.path.basename(input_path)
    filling_basename = os.path.basename(filling_path) if filling_path else None

    grid_parts = grid.split()
    if int(dimension) == 3 and len(grid_parts) == 2:
        grid_parts.append(grid_parts[-1])  # extend to a 3D grid, reusing the last value as z-split
    if int(dimension) == 3:
        computed_np = 1
        for g in grid_parts:
            computed_np *= int(g)
        np_procs = computed_np

    cmd_parts = [
        "mpirun", "--allow-run-as-root","--oversubscribe",
        "-np", str(np_procs),
        "./microsim_gp",
        input_basename,
    ]
    if filling_basename:
        cmd_parts.append(filling_basename)
    cmd_parts.append("DATA")
    cmd_parts.extend(grid_parts)
    if int(dimension) != 3:
        cmd_parts.append("0")
    print()
    print("=" * 60)
    print("RUNNING SIMULATION")
    print("=" * 60)
    print(f"  Command : {' '.join(cmd_parts)}")
    print(f"  Work dir: {work_dir}")
    print(f"  Output  : {data_dir}/")
    print("=" * 60)
    print()

    try:
        proc = subprocess.run(
            cmd_parts,
            cwd=work_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=3600,
        )
        retcode = proc.returncode
        print(proc.stdout.decode("utf-8", errors="replace"))
    except subprocess.TimeoutExpired:
        print("\nSimulation timed out after 1 hour.")
        retcode = -1
    except Exception as e:
        print(f"\nSimulation failed: {e}")
        retcode = -1

    if retcode == 0:
        data_files = os.listdir(data_dir) if os.path.isdir(data_dir) else []
        print()
        print("=" * 60)
        print("SIMULATION COMPLETE")
        print(f"  Output directory: {data_dir}/")
        print(f"  Files generated : {len(data_files)}")
        if data_files:
            for f in sorted(data_files)[:20]:
                print(f"    {f}")
            if len(data_files) > 20:
                print(f"    ... and {len(data_files) - 20} more")
        print("=" * 60)
    else:
        print(f"\nSimulation exited with code {retcode}")

    return retcode


# ============================================================================
# PIPELINE ORCHESTRATION
# ============================================================================

def run_stage(script_name, args, description):
    """Run one pipeline stage as a subprocess. Returns True on success."""
    cmd = [sys.executable, os.path.join(SCRIPT_DIR, script_name)] + args
    print()
    print("=" * 60)
    print(f"  STAGE: {description}")
    print("=" * 60)
    proc = subprocess.run(cmd)
    return proc.returncode == 0


def run_pipeline(prompt, base_dir, np_procs, grid, pipeline_dir, interactive):
    os.makedirs(pipeline_dir, exist_ok=True)
    spec_path = os.path.join(pipeline_dir, "spec.json")
    build_path = os.path.join(pipeline_dir, "build_result.json")
    validation_path = os.path.join(pipeline_dir, "validation_result.json")

    # --- Stage 1: translate ---
    stage1_args = ["--base-dir", base_dir, "-o", spec_path]
    if prompt:
        stage1_args = [prompt, "--base-dir", base_dir, "-o", spec_path]
    if not interactive:
        stage1_args.append("--non-interactive")
    if not run_stage("stage1_translator.py", stage1_args, "1/3 Translating request"):
        print("\n  Pipeline stopped: translation stage failed.")
        return 1

    # --- Stage 2: build ---
    stage2_args = [spec_path, "-o", build_path, "--base-dir", base_dir]
    if not interactive:
        stage2_args.append("--non-interactive")
    if not run_stage("stage2_builder.py", stage2_args, "2/3 Building simulation files"):
        print("\n  Pipeline stopped: build stage failed.")
        return 1

    # --- Stage 3: validate ---
    if not run_stage("stage3_validator.py",
                      [build_path, "-o", validation_path],
                      "3/3 Validating thermodynamic data"):
        print("\n  Pipeline stopped: validation failed. See issues above.")
        print(f"  (Re-run stage3_validator.py with --force to run anyway.)")
        return 1

    # --- Run the simulation ---
    with open(build_path) as f:
        build_result = json.load(f)

    return run_simulation(
        build_result["input_path"],
        build_result.get("filling_path"),
        np_procs=np_procs, grid=grid,
        work_dir=base_dir,
        system_name=build_result.get("system"),
        base_dir=base_dir,
        dimension=build_result.get("dimension", 2),
    )


# ============================================================================
# MAIN ENTRY POINT
# ============================================================================

def main():
    base_dir = os.environ.get("MICROSIM_DIR", "/app/microsim_gp")
    np_procs = int(os.environ.get("MPI_NP", "4"))
    grid = os.environ.get("MPI_GRID", "2 2")
    pipeline_dir = os.path.join(base_dir, ".pipeline")

    if "--help" in sys.argv or "-h" in sys.argv:
        print(__doc__)
        sys.exit(0)

    if len(sys.argv) > 1:
        prompt = " ".join(sys.argv[1:])
        sys.exit(run_pipeline(prompt, base_dir, np_procs, grid, pipeline_dir, interactive=False))
    else:
        # Interactive: stage1 itself prompts for input if no prompt is given.
        sys.exit(run_pipeline(None, base_dir, np_procs, grid, pipeline_dir, interactive=True))


if __name__ == "__main__":
    main()
