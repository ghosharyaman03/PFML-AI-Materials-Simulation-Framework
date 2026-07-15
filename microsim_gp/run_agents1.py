import autogen
import os
import re
import json
import copy
import sys

# ============================================================================
# 1. CONFIGURATION SCHEMA
# ============================================================================
#
# Any microsim_gp parameter can be passed as a dict entry. Keys are
# case-sensitive and must match the .in file syntax exactly.
#
# Value conventions:
#   scalar      → int | float | str          → KEY = value;
#   flat list   → [v1, v2, ...]              → KEY = {v1, v2, ...};
#   repeated    → [[...], [...], ...]        → one KEY = {...}; line per inner list
#
# Scientific notation that must survive round-tripping should be quoted:
#   "3e-8" stays as "3e-8" in the .in file, while 0.00000003 would not.

MULTI_LINE_KEYS = frozenset({
    "DIFFUSIVITY", "EIGEN_STRAIN", "VOIGT_ISOTROPIC",
    "BOUNDARY", "BOUNDARY_VALUE",
    "ceq", "cfill", "c_guess", "slopes", "A",
    "FILLCYLINDER", "FILLSPHERE", "FILLCUBE", "FILLRECT",
})

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
      "EIGEN_STRAIN", "VOIGT_ISOTROPIC"]),
    ("Boundary conditions",
     ["BOUNDARY", "BOUNDARY_VALUE"]),
    ("Model-specific parameters",
     ["ISOTHERMAL", "BINARY", "TERNARY", "DILUTE", "T",
      "WRITEFORMAT", "WRITEHDF5",
      "epsilon", "tau", "Tau"]),
    ("Anisotropy",
     ["Function_anisotropy", "Anisotropy_type", "dab",
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
]

DEFAULT_CONFIG = {
    "DIMENSION": 2,
    "MESH_X": 100, "MESH_Y": 100, "MESH_Z": 1,
    "DELTA_X": "1e-7", "DELTA_Y": "1e-7", "DELTA_Z": "1e-7", "DELTA_t": "1e-7",
    "NUMPHASES": 2, "NUMCOMPONENTS": 2,
    "NTIMESTEPS": 20000, "NSMOOTH": 10, "SAVET": 1000,
    "RESTART": 0, "STARTTIME": 0, "numworkers": 4, "TRACK_PROGRESS": 10,
    "COMPONENTS": ["Al", "Ni"],
    "PHASES": ["alpha", "liquid"],
    "GAMMA": [0.1], "R": 8.314, "V": "10e-6",
    "DIFFUSIVITY": [[1, 0, 0], [1, 1, "1e-9"]],
    "EIGEN_STRAIN": [
        [0, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0],
        [1, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0],
    ],
    "VOIGT_ISOTROPIC": [[0, 270, 187.5, 125.0], [1, 270, 187.5, 125.0]],
    "BOUNDARY": [
        ["phi", 1, 1, 1, 1, 0, 0], ["mu", 1, 1, 1, 1, 0, 0],
        ["c", 1, 1, 1, 1, 0, 0], ["T", 1, 1, 1, 1, 0, 0],
    ],
    "BOUNDARY_VALUE": [
        ["phi", 0, 0, 0, 0, 0, 0], ["mu", 0, 0, 0, 0, 0, 0],
        ["c", 0, 0, 0, 0, 0, 0], ["T", 0, 0, 0, 0, 0, 0],
    ],
    "ISOTHERMAL": 1, "BINARY": 1, "TERNARY": 0, "DILUTE": 0,
    "T": 1717, "WRITEFORMAT": "ASCII", "WRITEHDF5": 1,
    "epsilon": "4e-7", "tau": 1.31, "Tau": [0.28],
    "Function_anisotropy": 1, "Anisotropy_type": 4, "dab": [0.01],
    "Rotation_matrix": [0, 1, 0, 0, 45],
    "Function_W": 1, "Gamma_abc": [],
    "Shift": 0, "Shiftj": 30,
    "Writecomposition": 1,
    "Noise_phasefield": 1, "Amp_Noise_Phase": 0.1,
    "Equilibrium_temperature": 1718, "Filling_temperature": 1717,
    "Tempgrady": [856.95, 0.1, "10e-6", 0, "320e-6"],
    "Function_F": 4,
    "num_thermo_phases": 2, "tdbfname": "alzn_mey.tdb",
    "tdb_phases": ["FCC_A1", "LIQUID"], "phase_map": ["FCC_A1", "LIQUID"],
    "ceq": [
        [0, 0, 0.0944758], [0, 1, 0.102866],
        [1, 0, 0.102866], [1, 1, 0.102866],
    ],
    "cfill": [
        [0, 0, 0.0944758], [0, 1, 0.102866],
        [1, 0, 0.102866], [1, 1, 0.102866],
    ],
}

# Alloy presets — override dicts layered on top of DEFAULT_CONFIG.
# "_filling" is a reserved meta-key for filling-file parameters.
PRESETS = {
    "NiAlMo": {
        "NUMCOMPONENTS": 3,
        "DELTA_X": "3e-8", "DELTA_Y": "3e-8", "DELTA_Z": "3e-8",
        "DELTA_t": "220e-9", "NTIMESTEPS": 200000,
        "COMPONENTS": ["Al", "Mo", "Ni"],
        "BINARY": 0, "TERNARY": 1, "T": 1580,
        "DIFFUSIVITY": [[1, 0, 0, 0], [1, 1, "1e-9", "0.5e-9"]],
        "epsilon": "12e-8",
        "Equilibrium_temperature": 1590, "Filling_temperature": 1580,
        "ceq": [
            [0, 0, 0.100035, 0.0878765], [0, 1, 0.0998044, 0.168258],
            [1, 0, 0.0998044, 0.168258], [1, 1, 0.0998044, 0.168258],
        ],
        "cfill": [
            [0, 0, 0.100035, 0.0878765], [0, 1, 0.0998044, 0.168258],
            [1, 0, 0.0998044, 0.168258], [1, 1, 0.0998044, 0.168258],
        ],
        "_filling": {"FILLCYLINDER": [[0, 0, 0, 0, 0, 5]]},
    },
    "NiAlCr": {
        "NUMCOMPONENTS": 3,
        "MESH_X": 1000, "MESH_Y": 1000,
        "DELTA_X": "1.5e-8", "DELTA_Y": "1.5e-8", "DELTA_Z": "1.5e-8",
        "DELTA_t": "7.5e-8", "NTIMESTEPS": 2000000,
        "COMPONENTS": ["Al", "Mo", "Ni"],
        "BINARY": 0, "TERNARY": 1, "T": 1680,
        "DIFFUSIVITY": [[1, 0, 0.0, 0.0], [1, 1, "1.0e-9", "0.5e-9"]],
        "VOIGT_ISOTROPIC": [
            [0, "166.83e9", "115.5e9", "77e9"],
            [1, "166.83e9", "115.5e9", "77e9"],
        ],
        "epsilon": "6e-8", "Amp_Noise_Phase": 0.01,
        "Equilibrium_temperature": 1682, "Filling_temperature": 1680,
        "A": [[0, 1, 1, 1, 1], [1, 1, 1, 1, 1]],
        "ceq": [
            [0, 0, 0.100035, 0.0878765], [0, 1, 0.0998044, 0.168258],
            [1, 1, 0.0998044, 0.168258], [1, 0, 0.0998044, 0.168258],
        ],
        "cfill": [
            [0, 0, 0.100035, 0.0878765], [0, 1, 0.0998044, 0.168258],
            [1, 1, 0.0998044, 0.168258], [1, 0, 0.0998044, 0.168258],
        ],
        "c_guess": [
            [0, 0, 0.100098, 0.091893], [0, 1, 0.0990562, 0.177906],
            [1, 1, 0.0990562, 0.177906], [1, 0, 0.0990562, 0.177906],
        ],
    },
    "AlZn_dendrite": {
        "MESH_X": 200, "MESH_Y": 100, "MESH_Z": 250,
        "DELTA_X": "2e-8", "DELTA_Y": "2e-8", "DELTA_Z": "2e-8",
        "DELTA_t": "18e-9", "NTIMESTEPS": 50000, "SAVET": 10000,
        "TRACK_PROGRESS": 200, "numworkers": 512,
        "COMPONENTS": ["Al", "Zn"],
        "ISOTHERMAL": 0, "T": 857,
        "WRITEFORMAT": "BINARY",
        "ELASTICITY": 1,
        "BOUNDARY": [
            ["phi", 3, 3, 1, 1, 3, 3], ["mu", 3, 3, 1, 1, 3, 3],
            ["u", 3, 3, 1, 1, 3, 3], ["c", 3, 3, 1, 1, 3, 3],
            ["T", 3, 3, 1, 1, 3, 3], ["F", 3, 3, 1, 1, 3, 3],
        ],
        "BOUNDARY_VALUE": [
            ["phi", 0, 0, 0, 0, 0, 0], ["mu", 0, 0, 0, 0, 0, 0],
            ["c", 0, 0, 0, 0, 0, 0], ["T", 0, 0, 0, 0, 0, 0],
            ["F", 0, 0.2, 0, 0, 0, 0],
        ],
        "epsilon": "10e-8",
        "Rotation_matrix": [0, 1, 0, 0, 0],
        "Shift": 1, "Shiftj": 800,
        "Equilibrium_temperature": 870, "Filling_temperature": 857,
        "Tempgrady": [0.88, 0.12, 19200, 0, 0.010],
        "Function_F": 2,
        "A": [[0, 1], [1, 1]],
        "ceq": [
            [0, 0, 0.926], [0, 1, 0.817],
            [1, 1, 0.817], [1, 0, 0.817],
        ],
        "cfill": [
            [0, 0, 0.926], [0, 1, 0.926],
            [1, 1, 0.926], [1, 0, 0.926],
        ],
        "c_guess": [
            [0, 0, 0.92133], [0, 1, 0.81354],
            [1, 1, 0.81354], [1, 0, 0.81354],
        ],
        "slopes": [
            [0, 0, 772.34], [0, 1, 299.26],
            [1, 0, 299.26], [1, 1, 299.26],
        ],
    },
    "NiAl": {
        "NUMCOMPONENTS": 2,
        "COMPONENTS": ["Al", "Ni"],
        "BINARY": 1, "TERNARY": 0, "T": 1600,
        "NTIMESTEPS": 30000,
        "Equilibrium_temperature": 1610,
        "_filling": {"FILLRECT": [[20, 20, 0, 80, 80, 0, 1]]},
    },
    "NiNb": {
        "NUMCOMPONENTS": 2,
        "COMPONENTS": ["Nb", "Ni"],
        "BINARY": 1, "TERNARY": 0, "T": 1450,
        "NTIMESTEPS": 40000,
        "Equilibrium_temperature": 1460,
    },
}

# ============================================================================
# 2. TEMPLATE ENGINE — generate / parse .in files
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
    """Render a config dict into .in file content string."""
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
    """Render filling parameters into filling .in file content."""
    if not filling_config:
        return None
    lines = []
    for key, value in filling_config.items():
        if _is_repeated(key, value):
            for entry in value:
                lines.append(f"{key} = {_fmt_list(entry)};")
        elif isinstance(value, list):
            lines.append(f"{key} = {_fmt_list(value)};")
        else:
            lines.append(f"{key} = {_fmt(value)};")
    return "\n".join(lines) + "\n"


def parse_input_file(filepath):
    """Parse an existing .in file into a config dict."""
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
        return [_coerce(p.strip()) for p in inner.split(",")]
    return _coerce(raw)


def _coerce(s):
    try:
        if "." not in s and "e" not in s.lower():
            return int(s)
        return float(s)
    except (ValueError, AttributeError):
        return s

# ============================================================================
# 3. NATURAL LANGUAGE PROMPT PARSER
# ============================================================================

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
}


def parse_prompt(prompt):
    """Extract a preset name and numeric overrides from natural language.

    Returns (preset_name | None, overrides_dict, original_prompt).

    Examples:
        "Run NiAlMo at 1500K with 4 phases"
            → ("NiAlMo", {"T": 1500, "NUMPHASES": 4}, ...)
        "Set temperature to 1580 and timesteps to 100000"
            → (None, {"T": 1580, "NTIMESTEPS": 100000}, ...)
        "AlZn_dendrite with 200000 timesteps and dimension 3"
            → ("AlZn_dendrite", {"NTIMESTEPS": 200000, "DIMENSION": 3}, ...)
    """
    preset = None
    overrides = {}

    for name in sorted(PRESETS, key=len, reverse=True):
        if re.search(rf"\b{re.escape(name)}\b", prompt, re.IGNORECASE):
            preset = name
            break

    lookup = {}
    for alias, key in PARAM_ALIASES.items():
        lookup[alias] = key
    for key in DEFAULT_CONFIG:
        low = key.lower()
        if low not in lookup:
            lookup[low] = key

    num = r"(\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)"

    for alias in sorted(lookup, key=len, reverse=True):
        config_key = lookup[alias]
        if config_key in overrides:
            continue
        escaped = re.escape(alias)

        m = re.search(
            rf"\b{escaped}\b\s*(?:to|=|:|of)?\s*{num}\s*[kKcC]?\b",
            prompt, re.IGNORECASE,
        )
        if not m:
            m = re.search(
                rf"{num}\s*[kKcC]?\s+{escaped}\b",
                prompt, re.IGNORECASE,
            )
        if m:
            overrides[config_key] = _coerce(m.group(1))

    return preset, overrides, prompt

# ============================================================================
# 4. CONFIGURATION BUILDER
# ============================================================================

def build_config(base=None, preset=None, overrides=None, from_file=None):
    config = copy.deepcopy(DEFAULT_CONFIG)

    if base:
        config.update(copy.deepcopy(base))

    filling = None
    if preset and preset in PRESETS:
        preset_data = copy.deepcopy(PRESETS[preset])
        filling = preset_data.pop("_filling", None)
        config.update(preset_data)

    if from_file and os.path.exists(from_file):
        config.update(parse_input_file(from_file))

    if overrides:
        config.update(overrides)

    # --- ADD THIS SANITY PATCH TO AUTO-EXPAND ARRAYS ---
    num_phases = int(config.get("NUMPHASES", 2))
    current_phases = config.get("PHASES", [])
    
    # Correct the PHASES array size dynamically
    if len(current_phases) != num_phases:
        extended_phases = ["alpha", "liquid", "beta", "gamma", "delta", "epsilon", "zeta"]
        config["PHASES"] = extended_phases[:num_phases]
        while len(config["PHASES"]) < num_phases:
            config["PHASES"].append(f"phase_{len(config['PHASES'])}")

    return config, filling

# ============================================================================
# 5. FILE WRITER
# ============================================================================

def write_simulation_files(config, filling_config=None,
                           output_dir="/app/microsim_gp",
                           input_name="Input_generated.in",
                           filling_name="Filling_generated.in"):
    """Write .in files to disk and return their paths."""
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
# 6. AUTOGEN AGENTS
# ============================================================================

config_list = [
    {
        "model": "gemini-2.5-flash",
        "api_key": "AQ.Ab8RN6LRjgUt2x6pjzD3shG-yU04rjzNwCEeqNs5-W9bjUIssg",
        "base_url": "https://generativelanguage.googleapis.com/v1beta/openai/"
    }
]

ORCHESTRATOR_SYSTEM_MSG = """\
You are an automated expert engineer managing phase-field metallurgy simulations via `microsim_gp`.

The simulation controller has already generated one or more .in files for you.
Your job is to verify them, prepare the output directory, and execute the run.

Root path inside the container: /app/microsim_gp

WORKFLOW:
1. Cat the generated input file to verify it looks reasonable.
2. If a filling file was also generated, cat that too.
3. Create the DATA output directory if it does not exist:
       mkdir -p /app/microsim_gp/DATA
4. Run the simulation using the exact mpirun command provided to you.

MPI Execution Rules:
- With BOTH Input and Filling files:
    mpirun --allow-run-as-root -np <NP> /app/microsim_gp/microsim_gp <INPUT> <FILLING> ./DATA <GRID> 2>&1 | iconv -c -t UTF-8
- With ONLY an Input file:
    mpirun --allow-run-as-root -np <NP> /app/microsim_gp/microsim_gp <INPUT> ./DATA <GRID> 2>&1 | iconv -c -t UTF-8

CRITICAL: Always append `2>&1 | iconv -c -t UTF-8` to strip non-UTF-8 characters
that crash the AutoGen framework.

Once logs show the initialization phase passes and time steps start saving, output 'TERMINATE'.
"""

orchestrator = autogen.AssistantAgent(
    name="Orchestrator",
    system_message=ORCHESTRATOR_SYSTEM_MSG,
    llm_config={"config_list": config_list, "temperature": 0.1},
)

user_proxy = autogen.UserProxyAgent(
    name="UserProxy",
    human_input_mode="NEVER",
    max_consecutive_auto_reply=8,
    is_termination_msg=lambda x: x.get("content", "") and "TERMINATE" in x.get("content", ""),
    code_execution_config={
        "work_dir": "/app/microsim_gp",
        "use_docker": False,
    },
)

# ============================================================================
# 7. MAIN ENTRY POINT
# ============================================================================

def run(user_prompt=None, config_overrides=None, filling_overrides=None,
        base_config=None, from_file=None, preset=None,
        np_procs=4, grid="2 2",
        output_dir="/app/microsim_gp"):
    """Run a simulation with fully dynamic parameters.

    Three ways to supply parameters (all composable):

    1. Preset name:
         run(preset="NiAlMo")

    2. Dict overrides:
         run(config_overrides={"T": 1500, "NUMPHASES": 4})

    3. Natural language prompt (parsed automatically):
         run(user_prompt="Run NiAlMo at 1500K with 4 phases")

    Priority: defaults → preset → from_file → config_overrides → prompt overrides
    """
    detected_preset = None
    nl_overrides = {}
    if user_prompt:
        detected_preset, nl_overrides, _ = parse_prompt(user_prompt)

    active_preset = preset or detected_preset

    all_overrides = {**(config_overrides or {}), **nl_overrides}

    config, preset_filling = build_config(
        base=base_config, preset=active_preset,
        overrides=all_overrides, from_file=from_file,
    )

    filling = preset_filling or {}
    if filling_overrides:
        filling.update(filling_overrides)

    input_path, filling_path = write_simulation_files(
        config, filling or None, output_dir=output_dir,
    )

    print("=" * 60)
    print("SIMULATION CONFIGURATION")
    print("=" * 60)
    print(f"  Preset       : {active_preset or 'none (defaults)'}")
    if all_overrides:
        print(f"  Overrides    : {json.dumps(all_overrides, indent=15)}")
    if from_file:
        print(f"  Base file    : {from_file}")
    print(f"  Input file   : {input_path}")
    if filling_path:
        print(f"  Filling file : {filling_path}")
    print(f"  MPI procs    : {np_procs}")
    print(f"  MPI grid     : {grid}")
    print("=" * 60)

    if filling_path:
        mpi_cmd = (
            f"mpirun --allow-run-as-root -np {np_procs} "
            f"/app/microsim_gp/microsim_gp {input_path} {filling_path} "
            f"{grid} 2>&1 | iconv -c -t UTF-8"
        )
    else:
        mpi_cmd = (
            f"mpirun --allow-run-as-root -np {np_procs} "
            f"/app/microsim_gp/microsim_gp {input_path} "
            f"{grid} 2>&1 | iconv -c -t UTF-8"
        )

    # Note the indentation here (4 spaces) matching the rest of the run() function
    task_msg = (
        f"A simulation input file has been generated at: {input_path}\n"
        + (f"A filling file has been generated at: {filling_path}\n"
           if filling_path else "")
        + f"\nSteps:\n"
        f"1. Verify the input:  cat {input_path}\n"
        + (f"   Also verify:       cat {filling_path}\n" if filling_path else "")
        + f"2. Prepare output:    mkdir -p {output_dir}/DATA\n"
        f"3. Execute:           cd /app/microsim_gp && {mpi_cmd}\n\n"
        f"Monitor the execution until all time steps are finished. Once the simulation completes completely, output 'TERMINATE'." 
    )

    if user_prompt:
        task_msg = f"User request: \"{user_prompt}\"\n\n{task_msg}"

    print(f"\n--- Starting Agent Session ---\n")
    user_proxy.initiate_chat(orchestrator, message=task_msg)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        run(user_prompt=" ".join(sys.argv[1:]))
    else:
        run()
