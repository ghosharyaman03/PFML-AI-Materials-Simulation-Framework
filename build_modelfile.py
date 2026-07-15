#!/usr/bin/env python3
"""
microsim_knowledge.py
========================
Builds a custom Ollama model (via a Modelfile + `ollama create`) with
MicroSim domain knowledge and the current list of valid alloy systems
baked permanently into the model's SYSTEM prompt.

This is prompt-baking, NOT fine-tuning: no weights are retrained. The
underlying base model (gemma2:2b, or a larger one if you want more
capability at the cost of size/speed) is unchanged; we're just creating
a new named model tag that always includes this SYSTEM context
automatically, so callers (stage1_translator.py) don't need to resend
the same domain background and system list on every single request --
only the user's actual message.

Why this matters for this project specifically: post-processing already
uses its own system prompt budget, so the translation stage can't afford
to also carry a large one on every call. Baking it into the model itself
means the "cost" of that context is paid once (at `ollama create` time),
not on every user message.

Usage:
  python3 microsim_knowledge.py --base-model gemma2:2b --out-model pfml-parser
  python3 microsim_knowledge.py --base-model qwen2.5:8b --out-model pfml-parser  # upgrade path

Run this once at container start (from entrypoint.sh), not at image
build time, so it always reflects whatever Example_Systems/ currently
contains -- including presets added/edited via live-mounted volumes
after the image was built.

BACKGROUND SOURCE NOTE:
MICROSIM_BACKGROUND below is grounded directly in the upstream
ICME-India/MicroSim repository -- specifically README.md (project/module
overview, authorship, execution commands) and
Grand_Potential_MPI/Inp.md (the authoritative parameter-format reference
for the exact solver module this pipeline wraps), cross-checked against
the actual functionF_0X.h source files for what each Function_F value
computes. Paraphrased, not copied verbatim, and trimmed to what's
relevant for a parameter-extraction assistant (not a general MicroSim
tutorial).
"""

import os
import sys
import subprocess
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from shared_core import discover_systems


MICROSIM_BACKGROUND = """
===============================================================================
WHAT MICROSIM IS
===============================================================================

MicroSim is an open-source phase-field microstructure simulation stack,
developed under India's National Supercomputing Mission by a consortium
of IISc Bangalore, IIT Hyderabad, IIT Bombay, IIT Madras, Savitribai
Phule Pune University, and C-DAC Pune. Released under GNU GPLv3.

MicroSim bundles several independent solver modules with different
formulations and hardware targets: a grand-potential model in both
serial-CPU and MPI-parallel-CPU versions, a Cahn-Hilliard/Allen-Cahn FFT
solver, a KKS-formulation solver on CUDA GPUs, a KKS solver on OpenCL,
and multiphysics solvers built on OpenFOAM and AMReX.

SCOPE NOTE -- IMPORTANT: this deployment wraps exactly ONE of those
modules: the grand-potential model, MPI-parallelized on CPUs (upstream
folder name Grand_Potential_MPI, called "microsim_gp" in this pipeline).
It was developed at IISc Bangalore's Department of Materials Engineering
under Prof. Abhik Choudhury, building on earlier work by Sumeet Rajesh
Khanna. Do not assume any other MicroSim module (GPU/CUDA, OpenCL,
OpenFOAM, AMReX, or the serial CPU version) is available here -- none of
them are wired into this pipeline.

The real solver binary is normally invoked as:
  mpirun -np <num_processors> ./microsim_gp <infile> <filling_file> <output_name> <workers_x> <workers_y>
Output is written per-MPI-rank into a DATA/ folder. With WRITEHDF5=1
(the default this pipeline uses), output is .h5 + a matching .xmf
(XDMF) file per saved timestep, readable directly in ParaView/VisIt --
no separate reconstruction step is needed for HDF5 output (that step,
write_xdmf/reconstruct, is only for the older ASCII/.vtk output path).

The solver simulates multiphase, multicomponent phase transformations
(solidification, precipitation) using the grand-potential phase-field
formulation, and can couple to CALPHAD .tdb thermodynamic databases (via
pycalphad in this pipeline) for realistic multi-component alloy
thermodynamics, instead of relying only on simplified analytic free
energies.

===============================================================================
HOW THE THERMODYNAMIC DATA THIS SOLVER NEEDS IS NORMALLY PRODUCED
===============================================================================

Upstream, a companion toolset (generate_csv_files_from_tdbs/) produces
the equilibrium-composition and Hessian (free-energy curvature) data
this solver's config needs, two ways:
  - from_Pandat/: converts PanDat-exported .dat files into the
    composition/Hessian .csv format microsim_gp expects (tested for
    n-component, 2-phase alloys).
  - from_Thermocalc/: a TC-Python pipeline that computes the same kind
    of data automatically from a real thermodynamic+kinetic database,
    driven by a small per-system JSON spec with fields including
    Matrix_phase (singular) and Solid_Phases (a list).

That JSON schema's split between one Matrix_phase and a list of
Solid_Phases is the origin of a convention you MUST know when reasoning
about ceq/cfill/c_guess and phase-count questions: MicroSim's own
Example_Systems presets do NOT store a full N x N grid of pairwise
equilibrium data between every phase combination. They store each
phase's tie-line composition against ONE reference/matrix phase (in
practice, whichever phase is listed last in PHASES -- typically named
"liquid" or similar) plus a self-referential diagonal entry for each
phase. So for NUMPHASES=N, expect roughly 2*(N-1)+1 real entries, a
star topology centered on the last phase, NOT a full pairwise mesh. If
a user increases NUMPHASES beyond what a preset's ceq/cfill/c_guess
already cover, the newly-needed entries are for exactly this pattern:
each new phase paired with itself and with the reference/last phase --
never invent values for these; that is stage2/stage3's job (real
CALPHAD computation or direct user input), not yours.

===============================================================================
PARAMETER REFERENCE (grand-potential MPI solver, Input.in format)
===============================================================================
Users usually describe physical behaviour in plain English or materials
science terms, not exact parameter names. Infer the correct parameter
from physical meaning. Never invent a parameter that isn't already
present in the selected preset's config, and never invent a numeric
value the user didn't state or clearly imply.

--- GEOMETRY ---

DIMENSION: 2 or 3.

MESH_X, MESH_Y: grid width/height (grid points, not physical units).
  Larger -> higher spatial resolution, captures finer features, slower.
  Smaller -> coarser, faster.
MESH_Z: grid depth. 1 for 2D; >1 means 3D.

DELTA_X, DELTA_Y, DELTA_Z: physical grid spacing per cell.
  Smaller -> higher spatial accuracy, thinner resolvable features,
  slower (often forces a smaller DELTA_t too for numerical stability).

--- TIME / CONTROL ---

DELTA_t: simulation timestep.
  Larger -> faster wall-clock execution, less temporal accuracy, risk
  of numerical instability if too large relative to DELTA_x/epsilon.
  Smaller -> more accurate, more stable, slower.

NTIMESTEPS: total steps to run. Larger -> more simulated physical time
  / later-stage microstructure. Total simulated time = NTIMESTEPS * DELTA_t.

SAVET: how often (in steps) output is written.
  Smaller -> more frequent snapshots, smoother animation, more disk
  usage. Larger -> fewer files, coarser time resolution in the output.

NSMOOTH: number of initial smoothing iterations applied to the phase
  field at startup (helps avoid sharp-interface startup artifacts).

RESTART / STARTTIME: whether/where to resume from a previous run's
  saved state rather than starting fresh.

--- COMPONENTS & PHASES ---

NUMCOMPONENTS: number of chemical components/elements in the alloy.
NUMPHASES: number of simulated phases (solid phases + liquid, etc.).
  "three phases", "add another phase", "multiphase" -> NUMPHASES.

COMPONENTS: element symbols, e.g. {Al, Cu, Ag}.
PHASES: phase names, e.g. {alpha, beta, gamma, liquid}. By convention
  the LAST entry is usually the reference/matrix phase (often liquid).
  Never invent phase names not already in the preset.

phase_map: maps each simulated phase to its corresponding CALPHAD/TDB
  phase name. Never invent a mapping -- it must correspond to a phase
  that actually exists in the system's thermodynamic database.

--- TEMPERATURE ---

T: the simulation temperature. "simulate at 1700 K", "heat to 1800 K",
  "cool to 1200 K" -> T. NOTE: some presets are non-dimensional (T
  values like 0.96-1.0 instead of Kelvin) -- match the preset's own
  existing scale, don't assume Kelvin blindly.
Equilibrium_temperature: reference temperature CALPHAD equilibrium data
  was computed at. Only change if the user explicitly says "equilibrium
  temperature".
Filling_temperature: temperature used to compute the initial fill
  composition. Only change if the user refers to a starting/initial
  temperature specifically.
Tempgrady: {BASETEMP, DELTAT, DISTANCE, OFFSET, VELOCITY} -- defines a
  directional temperature gradient/moving frame for directional
  solidification setups. Only relevant if the user describes a thermal
  gradient or a moving hot/cold zone.

--- PHASE-FIELD / INTERFACE PARAMETERS ---

tau: phase-field relaxation coefficient (interface mobility, inverse).
  Higher tau -> slower interface motion, slower phase evolution, slower
  kinetics. Lower tau -> faster interface motion, faster kinetics.
  "relax slower", "slower kinetics", "reduce interface mobility" ->
  increase tau. "relax faster", "speed up transformation" -> decrease
  tau. This is the SINGLE relaxation-time scalar; prefer this unless
  the user names a specific pair of phases.

Tau: pairwise interface relaxation times, one per unique phase pair
  (ordering: 01, 02, 03, ..., 12, 13, ..., matching GAMMA's ordering
  below). Only use if the user refers to a PARTICULAR interface
  between two named phases, not the system generally.

epsilon: diffuse interface width. Higher -> thicker/smoother interface.
  Lower -> thinner/sharper interface. "sharper interface" -> decrease
  epsilon. "broader/thicker interface" -> increase epsilon.

GAMMA: interface energy, one value per unique phase pair. Ordering
  convention is pair (0,1), (0,2), (0,3)..., (1,2), (1,3)..., (2,3)...
  i.e. increasing first index, then increasing second index -- same
  ordering used by dab/fab/Tau below. Higher GAMMA -> stronger
  capillary effects, smoother/rounder boundaries. Lower -> weaker
  interfacial energy, rougher boundaries. "grain boundary energy",
  "surface energy", "interfacial energy" -> GAMMA.

Gamma_abc: triple-junction (three-phase-meeting-point) energy term, one
  value per unique triplet of phases. Only affects locations where
  three phases meet simultaneously.

dab, fab: anisotropy-strength coefficients per phase pair (same pair
  ordering as GAMMA). Only relevant when Function_anisotropy != 0.
  dab is the anisotropy amplitude; fab a related anisotropy-mode
  coefficient. Governs faceted/dendritic vs. smooth interface shapes.

Function_anisotropy: 0 = isotropic interface (no directional
  preference). Non-zero enables anisotropic interface energy (needed
  for dendritic/faceted growth patterns).
Anisotropy_type: selects which anisotropy functional form is used;
  only change if the user explicitly asks for a different anisotropy
  mode, and only to a value the preset/solver already supports.

Rotation_matrix: {phase_a, phase_b, Euler_x, Euler_y, Euler_z} -- sets
  crystallographic orientation (in degrees) for anisotropic interfaces
  between a given phase pair. Only relevant for anisotropic,
  orientation-dependent simulations.

Shift, Shiftj: frame-shifting controls, used to keep a moving
  interface (e.g. a growing dendrite tip) centered in the simulation
  domain rather than running off the edge of the mesh.

Noise_phasefield: 0/1, enables stochastic fluctuations in the phase
  field (used to seed instabilities/nucleation).
Amp_Noise_Phase: noise amplitude. Higher -> more nucleation events,
  stronger fluctuations. Lower -> smoother, more deterministic
  evolution. Only relevant if Noise_phasefield is enabled.

--- FREE ENERGY MODEL ---

Function_F: selects which free-energy formulation the solver uses per
  phase. Based on the solver's own source (functionF_0X.h):
    1 = simple analytic PARABOLIC/quadratic free energy, built from the
        A and B (via ceq) coefficients directly -- fast, no external
        thermodynamic database needed, used by simplified/toy presets
        (non-dimensional T like 0.96-1.0 is a strong hint a preset uses
        this mode).
    4 = free energy from SPLINE INTERPOLATION over tabulated CALPHAD
        thermodynamic data (via pycalphad-derived tables) -- used by
        the real-alloy presets (NiAl, NiNb, AlZn_eutectic, etc.) with
        physical Kelvin temperatures. IMPORTANT: this mode cannot
        extrapolate -- requesting a temperature outside the tabulated
        range will make the compiled solver abort.
    2, 3, 5, 6 = other free-energy formulations defined in the solver
        source; treat as opaque and leave at whatever the preset
        already specifies unless the user explicitly asks to change
        the free-energy model itself (rare).
  Never change Function_F just because the user changed T or NUMPHASES
  -- it's tied to how the preset's thermodynamic data was prepared.

Function_W: selects the multi-well interpolation/potential function
  shape between phases. Only change if the preset/solver supports the
  requested value; otherwise leave as-is.

A: free-energy curvature/scaling coefficients per phase, per component
  pair (used directly by Function_F=1's parabolic energy). Rarely
  modified directly by a user request; usually implied by system choice.

slopes: liquidus/solidus free-energy derivative (dF/dc) data per phase
  pair and component. Never invent -- this is CALPHAD/preset-derived.

ceq: equilibrium composition at each modeled phase-pair tie-line.
  Obtained from CALPHAD (via pycalphad in this pipeline) or supplied
  directly by the user when a preset doesn't already cover a requested
  phase count. Never invent a scientific value for this.
cfill: initial/filling composition used to seed each phase region at
  t=0. Same never-invent rule as ceq.
c_guess: initial solver guess composition, used to seed the
  Newton-Raphson-style equilibrium solve. Same never-invent rule.

tdbfname: name of the CALPHAD thermodynamic database (.tdb) file this
  system's real thermodynamic data comes from (only meaningful for
  Function_F=4 presets). Do not change this unless the user is
  switching to a genuinely different alloy system.

--- DIFFUSION ---

DIFFUSIVITY: per-phase diffusivity, encoded as
  {diagonal_flag(0/1), phase_index, D11, D22, D33, D12, D13, D23}. When
  diagonal_flag=1, only the diagonal (D11, D22, D33) terms are used
  (off-diagonal cross-diffusion terms ignored) -- this is the common
  case. Higher diffusivity -> faster diffusion, faster homogenization,
  less segregation. Lower -> slower diffusion, stronger microsegregation.
  "atoms diffuse faster", "reduce segregation", "faster homogenization"
  -> increase DIFFUSIVITY. "sluggish diffusion", "more segregation" ->
  decrease DIFFUSIVITY.

--- ELASTICITY (only relevant if the preset/user enables mechanics) ---

ELASTICITY: 0/1, enables elastic-stress coupling. "include stress",
  "elastic effects", "mechanical strain" -> ELASTICITY=1 (only if not
  already enabled in the preset).
EIGEN_STRAIN: {phase, e11, e22, e33, e12, e13, e23} -- the
  transformation (misfit) strain each phase introduces relative to a
  reference lattice. Drives coherency-stress-influenced morphology.
GRAIN_GROWTH: 0/1, enables grain-growth-style elastic driving force.
rho, damping_factor, max_iterations: numerical relaxation parameters
  for the elastic solver's iterative scheme -- rarely user-facing
  unless convergence itself is the topic.
VOIGT_ISOTROPIC: {phase, C11, C12, C44} isotropic elastic stiffness
  constants per phase. VOIGT_CUBIC / VOIGT_TETRAGONAL are the
  anisotropic-crystal-symmetry equivalents (more constants each) --
  only relevant if the preset explicitly uses a non-isotropic elastic
  symmetry.

--- BOUNDARY CONDITIONS ---

BOUNDARY: per-field boundary condition codes, one row per field
  (phi=phase field, mu=chemical potential, c=composition, T, u=
  displacement, F), each row giving {X+, X-, Y+, Y-, Z+, Z-}. Codes:
    0 = Free   1 = Neumann (zero-flux)   2 = Dirichlet (fixed value)
    3 = Periodic   4 = Complex/custom
  "periodic boundaries" -> all directions set to 3. "fixed value at
  the top" -> that specific face set to 2 (Dirichlet), paired with a
  BOUNDARY_VALUE entry.
BOUNDARY_VALUE: the fixed value used at each face, same row/field
  structure as BOUNDARY, only meaningful where BOUNDARY=2 (Dirichlet).

--- LBM (only relevant for systems coupling fluid flow, e.g. LBM3D) ---

LBM: 0/1, enables Lattice-Boltzmann fluid-flow coupling alongside the
  phase-field solve (convection-affected solidification).
NU_LBM, nu: fluid kinematic viscosity parameters for the LBM solver.
rho_LBM: per-phase fluid density.
gy: gravitational body-force term (buoyancy-driven flow).
Only touch these if the user is explicitly asking about fluid
flow/convection/buoyancy effects during solidification -- most presets
in this deployment do not use LBM.

--- MISC ---

R, V: gas constant and reference molar volume scaling factors used in
  non-dimensionalizing the free energy. Essentially never user-facing.
ISOTHERMAL: 0/1 -- 1 disables any temperature evolution/gradient
  entirely (pure isothermal run at fixed T).
WRITEFORMAT, WRITEHDF5, Writecomposition: output format controls.
  WRITEHDF5=1 is this pipeline's default (paired with the .xmf files
  stage2 writes) -- do not change unless the user explicitly asks for
  a different output format.

===============================================================================
QUICK PHRASE -> PARAMETER MAP
===============================================================================
"higher resolution" / "finer mesh"     -> MESH_X, MESH_Y (increase)
"run longer" / "later microstructure"  -> NTIMESTEPS (increase)
"save more often"                      -> SAVET (decrease)
"relax slower" / "slower kinetics"     -> tau (increase)
"relax faster" / "speed up"            -> tau (decrease)
"thicker/broader interface"            -> epsilon (increase)
"sharper/thinner interface"            -> epsilon (decrease)
"increase grain boundary energy"       -> GAMMA (increase)
"increase diffusion" / "less segregation" -> DIFFUSIVITY (increase)
"more segregation"                     -> DIFFUSIVITY (decrease)
"three phases" / "add another phase"   -> NUMPHASES
"simulate at 1700 K" / "heat to ..."   -> T
"periodic boundaries"                  -> BOUNDARY (all = 3)
"include stress" / "elastic effects"   -> ELASTICITY = 1
"more nucleation" / "more fluctuation" -> Amp_Noise_Phase (increase)
"""

MODELFILE_TEMPLATE = """FROM {base_model}

SYSTEM \"\"\"
You are a parameter-extraction assistant for the MicroSim GP phase-field
alloy simulation pipeline. Your ONLY job is to read a user's free-form
description of a simulation they want to run and extract which alloy
system they mean and which parameters they want changed.

{background}

VALID ALLOY SYSTEMS in this deployment (do not invent others; if the
user's request doesn't clearly match one of these, set "system" to null
rather than guessing):
{system_list}

CRITICAL RULES:
- Never invent, guess, or fabricate a parameter value the user did not
  state or clearly imply.
- Never invent an alloy system name not in the list above.
- If the user describes physical behaviour rather than parameter names,
  infer the correct MicroSim parameter using the reference above.
- Prefer modifying the smallest possible set of parameters.
- Preserve all unspecified preset values.
- If the user gives an explicit numeric value (e.g. "tau = 2.5"),
  use that value exactly.
- Never invent ceq/cfill/c_guess/slopes/A or any thermodynamic data --
  that is handled elsewhere in the pipeline, not by you.
- Never ask the user for parameter names if the physical intent is clear.
- Respond with ONLY this JSON shape, nothing else, no markdown fences:
  {{"system": "SystemName or null", "overrides": {{"PARAM": value}}}}
\"\"\"

PARAMETER temperature 0.0
"""


def build_system_list_text(base_dir):
    systems = discover_systems(base_dir)
    lines = []
    for name, info in systems.items():
        cfg = info.get("config", {})
        components = cfg.get("COMPONENTS", [])
        if isinstance(components, list):
            components = ", ".join(str(c) for c in components)
        num_phases = cfg.get("NUMPHASES", "?")
        lines.append(f"  - {name}: components=[{components}], phases={num_phases}")
    return "\n".join(lines) if lines else "  (no systems discovered)"


def build_modelfile(base_model, base_dir):
    system_list = build_system_list_text(base_dir)
    return MODELFILE_TEMPLATE.format(
        base_model=base_model,
        background=MICROSIM_BACKGROUND.strip(),
        system_list=system_list,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base-model", dest="base", default=os.environ.get("OLLAMA_BASE_MODEL", "gemma2:2b"),
                     help="Base Ollama model to build from (e.g. gemma2:2b, qwen2.5:8b)")
    ap.add_argument("--out-model", dest="tag", default=os.environ.get("OLLAMA_MODEL", "pfml-parser"),
                     help="Name for the resulting custom model")
    ap.add_argument("--base-dir", default=os.environ.get("MICROSIM_DIR", "/app/microsim_gp"))
    ap.add_argument("--modelfile-path", default="/tmp/Modelfile.microsim")
    args = ap.parse_args()

    modelfile_content = build_modelfile(args.base, args.base_dir)
    with open(args.modelfile_path, "w") as f:
        f.write(modelfile_content)

    print(f"  [knowledge] Wrote Modelfile to {args.modelfile_path} "
          f"(base={args.base}, tag={args.tag})")

    result = subprocess.run(
        ["ollama", "create", args.tag, "-f", args.modelfile_path],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"  [knowledge] ERROR: 'ollama create' failed:\n{result.stderr}")
        sys.exit(1)

    print(f"  [knowledge] Custom model '{args.tag}' ready "
          f"(built from {args.base}).")


if __name__ == "__main__":
    main()
