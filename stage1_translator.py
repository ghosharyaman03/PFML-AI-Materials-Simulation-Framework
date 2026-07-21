#!/usr/bin/env python3
"""
stage1_translator.py
=====================
Turns a user's free-form paragraph into a structured simulation spec.

Supports two modes:
  1. Existing system: matches to a preset in Example_Systems/
  2. New system: creates a custom config from a template

For new systems, asks the user for basic parameters (components, phases,
mesh, timesteps, temperature) and auto-pads all thermodynamic values
from a template system.

Usage:
  Interactive:  python3 stage1_translator.py
  One-shot:     python3 stage1_translator.py "AlZn_eutectic with 6 phases..." -o spec.json

Env vars:
  MICROSIM_DIR    Path to microsim_gp directory (default: /app/microsim_gp)
  OLLAMA_HOST     Base URL of the Ollama server (default: http://localhost:11434)
  OLLAMA_MODEL    Model to use for parsing (default: pfml-parser)
"""

import os
import re
import sys
import json
import argparse
import urllib.request
import urllib.error

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from shared_core import (
    SECTIONS, PARAM_ALIASES, discover_systems, show_systems, _coerce,
)

OLLAMA_HOST = os.environ.get("OLLAMA_HOST") or "http://localhost:11434"
OLLAMA_MODEL = os.environ.get("OLLAMA_MODEL", "pfml-parser")

MAX_IDENTIFICATION_ATTEMPTS = 6

FIELDS_WORTH_ASKING = [
    ("system", "Which alloy system do you want to simulate?"),
    ("NUMPHASES", "How many phases should the simulation have?"),
    ("MESH_X", "What grid/mesh width (MESH_X) do you want?"),
    ("MESH_Y", "What grid/mesh height (MESH_Y) do you want?"),
    ("NTIMESTEPS", "How many timesteps should it run for?"),
]

NEW_SYSTEM_QUESTIONS = [
    ("COMPONENTS", "What are the chemical components/elements? (comma-separated, e.g. Al, Co)"),
    ("NUMPHASES", "How many phases?"),
    ("MESH_X", "Grid width (MESH_X)?"),
    ("MESH_Y", "Grid height (MESH_Y)?"),
    ("NTIMESTEPS", "How many timesteps?"),
    ("T", "Simulation temperature (K)?"),
    ("DIMENSION", "Dimension (2 or 3)?"),
]

ELEMENT_KEYWORDS = {
    "aluminum": "Al", "aluminium": "Al",
    "zinc": "Zn",
    "nickel": "Ni",
    "copper": "Cu",
    "molybdenum": "Mo",
    "niobium": "Nb",
    "chromium": "Cr",
    "cobalt": "Co",
    "iron": "Fe",
    "titanium": "Ti",
    "tungsten": "W",
    "tantalum": "Ta",
}

# All recognized element symbols (uppercase and lowercase)
ALL_ELEMENT_SYMBOLS = {
    "AL", "ZN", "NI", "CU", "MO", "NB", "CR", "CO",
    "FE", "TI", "W", "TA",
}

SYSTEM_NAME_ALIASES = {
    "NiAlCr": ["chromium", "chromate", "cr alloy"],
}
NUMBER_WORDS = {
    "one": 1, "two": 2, "three": 3, "four": 4, "five": 5,
    "six": 6, "seven": 7, "eight": 8, "nine": 9, "ten": 10,
}
GENERIC_NAME_TOKENS = {"model", "example", "examples", "system", "alloy"}

UNCHANGED_SIGNAL_WORDS = {
    "NUMPHASES": ["phase", "phases"],
    "MESH_X": ["mesh"],
    "MESH_Y": ["mesh"],
    "NTIMESTEPS": ["timestep", "timesteps", "step", "steps"],
}
UNCHANGED_PATTERN = r"(same|unchanged|as[- ]is|no change|keep it)"
QUIT_COMMANDS = {"quit", "exit", "q", "bye", "goodbye"}
MAX_GARBAGE_ATTEMPTS = 3

SIMULATION_DOMAIN_TERMS = [
    "mesh", "phase", "timestep", "temperature", "grid",
    "component", "dimension", "delta", "epsilon", "gamma",
    "diffusivity", "elasticity", "boundary", "noise",
    "anisotropy", "filling", "equilibrium", "save", "restart",
    "molar", "volume", "damping", "stress", "concentration",
    "composition", "simulation", "simulate", "alloy",
    "solidification", "precipitation", "microstructure",
    "eutectic", "peritectic", "spinodal", "nucleation",
    "binary", "ternary", "isothermal", "grain", "domain",
    "relaxation", "interface", "energy", "free energy",
    "function", "parabolic", "spline", "calphad", "thermodynamic",
    "database", "tdb", "matrix", "liquid", "solid", "fcc", "bcc",
    "hcp", "cubic", "input", "filling", "run", "execute",
]


def classify_input(prompt):
    """Fast regex-based filter to distinguish simulation requests from garbage.
    Returns 'valid' if the input could be a simulation request, 'garbage' otherwise.
    Handles full names ('aluminum'), symbols ('al'), and CamelCase ('NiAl').
    Rejects garbage values like '2000dx' or 'banana'."""
    p = prompt.strip().lower()
    if len(p) < 2:
        return "garbage"

    # 1. Element full names (aluminum, cobalt, nickel, etc.)
    for word in ELEMENT_KEYWORDS:
        if re.search(rf"\b{re.escape(word)}\b", p):
            return "valid"

    # 2. Element symbols: 1-3 letter words that are known symbols
    #    ('al', 'co', 'ni' → valid; 'aa', 'zz' → garbage)
    words = re.findall(r"\b[a-zA-Z]{1,3}\b", p)
    for w in words:
        if w.upper() in ALL_ELEMENT_SYMBOLS:
            return "valid"

    # 3. CamelCase system names like NiAl, AlZn, NiAlCr
    if re.search(r"[A-Z][a-z]+[A-Z]", prompt):
        return "valid"

    # 4. Underscore names like AlZn_eutectic
    if re.search(r"\w+_\w+", p):
        return "valid"

    # 5. Parameter keyword followed by = or : and a CLEAN number
    #    'phases=3' → valid, 'mesh=2000dx' → garbage, 'mesh_x=200' → valid
    if re.search(r"(?:mesh|phase|timestep|temperature|dimension|delta|epsilon|gamma|num)\w*\s*[=:]\s*\d+(?!\w)", p):
        return "valid"

    # 6. Parameter keyword followed by a CLEAN number
    #    'phases 3' → valid, 'mesh 2000' → valid, 'mesh 2000dx' → garbage
    if re.search(r"(?:mesh|phase|timestep|temperature|dimension|num)\w*\s+\d+(?!\w)", p):
        return "valid"

    # 7. Temperature with explicit unit: '1500K', '1500 k' (NOT bare '200')
    if re.search(r"\d+\s*[kK]\b", p) and re.search(r"(?:temp|t\s*=|temperature)", p):
        return "valid"

    return "garbage"


def check_input_quality(prompt):
    """Use Ollama to generate a friendly rejection message for garbage input.
    Always returns (False, message) -- this is only called after classify_input
    already determined the input is garbage, so we just need the AI to explain
    to the user what went wrong."""
    raw = ask_ollama(
        f"A user is interacting with a materials science simulation assistant. "
        f"The assistant helps users set up alloy phase-field simulations.\n\n"
        f"The user typed: '{prompt}'\n\n"
        f"This input was NOT recognized as a valid simulation request. "
        f"Generate a short, friendly message (1-2 sentences) telling the user "
        f"their input was not understood and asking them to describe the alloy "
        f"system or simulation they want to run. For example they could type "
        f"'aluminum cobalt' or 'NiAl with 6 phases'. "
        f'Respond with ONLY this JSON: {{"msg": "<your message>"}}'
    )
    if raw:
        try:
            text = re.sub(r"```json\s*|```\s*", "", raw).strip()
            parsed = json.loads(text)
            return False, parsed.get("msg", "")
        except Exception:
            pass
    return False, ""


def detect_unchanged_fields(prompt):
    skipped = set()
    for key, topic_words in UNCHANGED_SIGNAL_WORDS.items():
        for word in topic_words:
            if re.search(rf"\b{word}\b[^.]{{0,25}}\b{UNCHANGED_PATTERN}\b", prompt, re.IGNORECASE) or \
               re.search(rf"\b{UNCHANGED_PATTERN}\b[^.]{{0,25}}\b{word}\b", prompt, re.IGNORECASE):
                skipped.add(key)
                break
    return skipped


def extract_components_from_prompt(prompt):
    """Extract element symbols from the prompt.
    Matches both full names ('aluminum' → 'Al') and symbols ('al' → 'Al')."""
    components = []
    # Full names (aluminum, cobalt, nickel, etc.)
    for word, sym in ELEMENT_KEYWORDS.items():
        if re.search(rf"\b{re.escape(word)}\b", prompt, re.IGNORECASE):
            if sym not in components:
                components.append(sym)
    # Short symbols (al, co, ni, etc.) — only 1-2 letter words to avoid
    # matching random English words like 'in', 'on', 'at'
    # Map uppercase back to mixed-case from ELEMENT_KEYWORDS
    _sym_to_mixed = {v.upper(): v for v in ELEMENT_KEYWORDS.values()}
    words = re.findall(r"\b([A-Za-z]{1,2})\b", prompt)
    for w in words:
        sym = w.upper()
        if sym in ALL_ELEMENT_SYMBOLS:
            mixed = _sym_to_mixed[sym]
            if mixed not in components:
                components.append(mixed)
    return components


# ============================================================================
# InputCollector — validated per-field user input
# ============================================================================

VALID_ELEMENT_SYMBOLS = {v.upper() for v in ELEMENT_KEYWORDS.values()}


class InputCollector:
    """Handles all interactive user input with per-field validation.
    Every field that asks the user for a value validates it before accepting."""

    @staticmethod
    def _ask(prompt_text):
        """Raw input with quit/EOF handling. Returns (text, should_quit)."""
        return _get_user_input(prompt_text)

    @staticmethod
    def get_number(prompt, min_val=None, max_val=None, allow_float=True, default=None):
        """Ask for a numeric value. Reprompts until a valid number is entered.
        Returns the number, or default if user skips, or None if user quits."""
        while True:
            text, quit = InputCollector._ask(prompt)
            if quit:
                return None
            if not text:
                if default is not None:
                    return default
                print("  This field requires a value. Please enter a number.")
                continue
            try:
                val = float(text)
                if not allow_float:
                    if val != int(val):
                        print(f"  '{text}' is not a whole number. Please enter an integer (no decimals).")
                        continue
                    val = int(val)
                if min_val is not None and float(val) < min_val:
                    print(f"  Value must be at least {min_val}. Try again.")
                    continue
                if max_val is not None and float(val) > max_val:
                    print(f"  Value must be at most {max_val}. Try again.")
                    continue
                return val
            except ValueError:
                print(f"  '{text}' is not a valid number. Please enter a numeric value.")

    @staticmethod
    def get_int(prompt, min_val=None, max_val=None, default=None):
        """Ask for an integer value. Reprompts until valid."""
        while True:
            result = InputCollector.get_number(
                prompt, min_val=min_val, max_val=max_val,
                allow_float=False, default=default,
            )
            if result is None:
                return None
            return int(result)

    @staticmethod
    def get_components(prompt):
        """Ask for alloy components. Validates each is a known element symbol.
        Returns list of element symbols, or None if user quits."""
        while True:
            text, quit = InputCollector._ask(prompt)
            if quit:
                return None
            if not text:
                print("  Please enter at least one component (e.g. Al, Co).")
                continue
            parts = [c.strip().upper() for c in text.split(",") if c.strip()]
            if not parts:
                print("  Please enter at least one component (e.g. Al, Co).")
                continue
            invalid = [c for c in parts if c not in VALID_ELEMENT_SYMBOLS]
            if invalid:
                print(f"  Unknown element(s): {', '.join(invalid)}.")
                print(f"  Use valid symbols like: {', '.join(sorted(VALID_ELEMENT_SYMBOLS))}")
                continue
            return parts

    @staticmethod
    def get_yes_no(prompt, default=None):
        """Ask a yes/no question. Returns True/False, or default if skipped."""
        while True:
            text, quit = InputCollector._ask(prompt)
            if quit:
                return None
            text_lower = text.lower().strip()
            if text_lower in ("y", "yes"):
                return True
            if text_lower in ("n", "no"):
                return False
            if not text and default is not None:
                return default
            print("  Please enter 'y' or 'n'.")


def collect_new_system_info(interactive=True):
    """Ask user for basic parameters to define a new alloy system.
    Uses InputCollector for per-field validation — every numeric field
    reprompts until a valid number is entered.
    Returns None if user quits at any point."""
    if not interactive:
        return None

    info = {}

    components = InputCollector.get_components(
        "  What are the chemical components/elements? (comma-separated, e.g. Al, Co): "
    )
    if components is None:
        return None
    info["COMPONENTS"] = components
    info["NUMCOMPONENTS"] = len(components)

    n_comp = len(components)
    default_phases = n_comp + 1

    num_phases = InputCollector.get_int(
        f"  How many phases? (default {default_phases}): ",
        min_val=1, max_val=50, default=default_phases,
    )
    if num_phases is None:
        return None
    info["NUMPHASES"] = int(num_phases)

    mesh_x = InputCollector.get_int(
        "  Grid width (MESH_X)? (default 100): ",
        min_val=1, max_val=10000, default=100,
    )
    if mesh_x is None:
        return None
    info["MESH_X"] = int(mesh_x)

    mesh_y = InputCollector.get_int(
        "  Grid height (MESH_Y)? (default 100): ",
        min_val=1, max_val=10000, default=100,
    )
    if mesh_y is None:
        return None
    info["MESH_Y"] = int(mesh_y)

    timesteps = InputCollector.get_int(
        "  How many timesteps? (default 100000): ",
        min_val=1, max_val=100000000, default=100000,
    )
    if timesteps is None:
        return None
    info["NTIMESTEPS"] = int(timesteps)

    temp = InputCollector.get_number(
        "  Simulation temperature in K? (default 1000): ",
        min_val=1, max_val=100000, default=1000,
    )
    if temp is None:
        return None
    info["T"] = temp

    dim = InputCollector.get_int(
        "  Dimension (2 or 3)? (default 2): ",
        min_val=2, max_val=3, default=2,
    )
    if dim is None:
        return None
    info["DIMENSION"] = int(dim)

    return info


# ============================================================================
# Regex fallback parser
# ============================================================================

def parse_prompt_regex(prompt, system_names=None, systems=None):
    system_names = system_names or []
    preset = None
    overrides = {}

    for name in sorted(system_names, key=len, reverse=True):
        pattern = re.escape(name).replace("_", r"[_ ]?")
        if re.search(rf"\b{pattern}\b", prompt, re.IGNORECASE):
            preset = name
            break

    if not preset:
        for name, aliases in SYSTEM_NAME_ALIASES.items():
            if name in system_names and any(
                re.search(rf"\b{re.escape(alias)}\b", prompt, re.IGNORECASE) for alias in aliases
            ):
                preset = name
                break

    if not preset and systems:
        mentioned = {
            sym.upper() for word, sym in ELEMENT_KEYWORDS.items()
            if re.search(rf"\b{word}\b", prompt, re.IGNORECASE)
        }
        if mentioned:
            candidates = []
            for name in system_names:
                comps = systems[name]["config"].get("COMPONENTS", [])
                if isinstance(comps, str):
                    comps = [comps]
                comp_set = {str(c).upper() for c in comps}
                if mentioned == comp_set:
                    candidates.append(name)

            if len(candidates) == 1:
                preset = candidates[0]
            elif len(candidates) > 1:
                for name in candidates:
                    tokens = [
                        t for t in re.split(r"[_ ]+", name)
                        if len(t) > 3 and t.lower() not in GENERIC_NAME_TOKENS
                    ]
                    if any(re.search(rf"\b{re.escape(t)}\b", prompt, re.IGNORECASE) for t in tokens):
                        preset = name
                        break

    mesh_xy_match = re.search(r"\bmesh\b.{0,15}?(\d+)\s*(?:[x\u00d7]|by)\s*(\d+)", prompt, re.IGNORECASE)
    if mesh_xy_match:
        overrides["MESH_X"] = _coerce(mesh_xy_match.group(1))
        overrides["MESH_Y"] = _coerce(mesh_xy_match.group(2))
    if "NUMPHASES" not in overrides:
        words_pattern = "|".join(NUMBER_WORDS.keys())
        phase_word_match = re.search(rf"\b({words_pattern})\b\s+phases?\b", prompt, re.IGNORECASE)
        if phase_word_match:
            overrides["NUMPHASES"] = NUMBER_WORDS[phase_word_match.group(1).lower()]
    temp_k_match = re.search(r"\b(\d+(?:\.\d+)?)\s*[kK]\b", prompt)
    if temp_k_match:
        overrides["T"] = _coerce(temp_k_match.group(1))

    lookup = dict(PARAM_ALIASES)
    all_known_keys = set()
    for _, keys in SECTIONS:
        all_known_keys.update(keys)
    for key in sorted(all_known_keys):
        lookup.setdefault(key.lower(), key)

    single_num = r"\d+(?:\.\d+)?(?:[eE][+-]?\d+)?"
    num_list = rf"({single_num}(?:\s*,\s*{single_num})*)"

    for alias in sorted(lookup, key=len, reverse=True):
        config_key = lookup[alias]
        if config_key in overrides:
            continue
        escaped = re.escape(alias)
        m = re.search(
            rf"\b{escaped}\b\s*(?:to|=|:|of|values?)?\s*{num_list}\s*[kKcC]?\b",
            prompt, re.IGNORECASE,
        ) or re.search(
            rf"{num_list}\s*[kKcC]?\s+{escaped}\b", prompt, re.IGNORECASE,
        )
        if m:
            raw = m.group(1)
            overrides[config_key] = (
                [_coerce(p.strip()) for p in raw.split(",")] if "," in raw
                else _coerce(raw)
            )

    return preset, overrides


# ============================================================================
# Ollama-assisted parser
# ============================================================================

def _ollama_available():
    try:
        req = urllib.request.Request(f"{OLLAMA_HOST}/api/tags")
        with urllib.request.urlopen(req, timeout=2) as resp:
            return resp.status == 200
    except Exception:
        return False


def ask_ollama(user_message, system_prompt=None):
    payload = {
        "model": OLLAMA_MODEL,
        "prompt": user_message,
        "stream": False,
        "format": "json",
        "options": {"temperature": 0.0},
    }
    if system_prompt:
        payload["system"] = system_prompt
    try:
        req = urllib.request.Request(
            f"{OLLAMA_HOST}/api/generate",
            data=json.dumps(payload).encode("utf-8"),
            headers={"Content-Type": "application/json"},
        )
        with urllib.request.urlopen(req, timeout=60) as resp:
            body = json.loads(resp.read().decode("utf-8"))
            return body.get("response", "")
    except (urllib.error.URLError, TimeoutError, Exception):
        return None


def _extract_numbers(text):
    matches = re.findall(r"-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?", text)
    nums = set()
    for m in matches:
        try:
            nums.add(float(m))
        except ValueError:
            pass
    return nums


def _is_grounded(value, number_pool, tol=1e-9):
    if isinstance(value, bool):
        return True
    if isinstance(value, (int, float)):
        return any(abs(float(value) - n) < tol for n in number_pool)
    if isinstance(value, list):
        return all(_is_grounded(v, number_pool, tol) for v in value)
    return True


def parse_prompt_llm(prompt, system_names, known_overrides=None):
    if not _ollama_available():
        return None, {}, []

    all_params = set()
    for _, keys in SECTIONS:
        all_params.update(keys)
    all_params.update(PARAM_ALIASES.values())

    raw = ask_ollama(prompt)
    if not raw:
        return None, {}, []

    try:
        text = re.sub(r"```json\s*|```\s*", "", raw).strip()
        parsed = json.loads(text)
    except Exception:
        return None, {}, []

    system_lookup = {s.lower(): s for s in system_names}
    raw_system = parsed.get("system")
    resolved_system = (
        system_lookup.get(str(raw_system).strip().lower()) if raw_system else None
    )

    coerced = {}
    for k, v in (parsed.get("overrides") or {}).items():
        if k not in all_params:
            continue
        coerced[k] = _coerce(v) if isinstance(v, str) else v

    number_pool = _extract_numbers(prompt)
    grounded = {}
    mentioned_no_value = []
    for k, v in coerced.items():
        if _is_grounded(v, number_pool):
            grounded[k] = v
        else:
            print(f"  [translator] Discarding {k}={v} suggested by the model -- "
                  f"no matching number found in what you typed (looks fabricated).")
            mentioned_no_value.append(k)

    return resolved_system, grounded, mentioned_no_value


# ============================================================================
# Core
# ============================================================================

def merge_overrides(base, new):
    merged = dict(base)
    for k, v in new.items():
        merged.setdefault(k, v)
    return merged


def missing_fields(system, overrides, skipped=None):
    skipped = skipped or set()
    missing = []
    for key, question in FIELDS_WORTH_ASKING:
        if key == "system":
            if not system:
                missing.append((key, question))
        elif key not in overrides and key not in skipped:
            missing.append((key, question))
    return missing


def _get_user_input(prompt_text=">>> "):
    """Get input from user, handling EOF/KeyboardInterrupt and quit commands.
    Returns (text, should_quit)."""
    try:
        text = input(prompt_text).strip()
    except (EOFError, KeyboardInterrupt):
        print("\n  Goodbye.")
        return "", True
    if text.lower() in QUIT_COMMANDS:
        return text, True
    return text, False


def _try_parse(prompt, system_names, systems):
    """Parse a prompt through regex + LLM. Returns (system, overrides)."""
    regex_system, regex_overrides = parse_prompt_regex(prompt, system_names, systems=systems)
    llm_system, llm_overrides, _ = parse_prompt_llm(prompt, system_names)
    system = regex_system or llm_system
    overrides = merge_overrides(regex_overrides, llm_overrides)
    return system, overrides


def _offer_new_system(system_names):
    """Ask user if they want to create a new system. Returns spec dict or None."""
    print(f"\n  Available systems: {', '.join(system_names)}")
    text, quit = _get_user_input(
        "  Would you like to create a NEW alloy system? [y/N] "
    )
    if quit:
        return None

    if text.lower() in ("y", "yes"):
        new_info = collect_new_system_info(True)
        if new_info is None:
            return None
        if new_info.get("COMPONENTS"):
            n_comp = new_info.get("NUMCOMPONENTS", 2)
            new_info.setdefault("NUMPHASES", n_comp + 1)
            new_info.setdefault("MESH_X", 100)
            new_info.setdefault("MESH_Y", 100)
            new_info.setdefault("NTIMESTEPS", 100000)
            new_info.setdefault("T", 1000)
            new_info.setdefault("DIMENSION", 2)
            return {"system": "NEW", "overrides": new_info, "raw_prompt": ""}
    return None


NUMERIC_FIELDS = {"NUMPHASES", "MESH_X", "MESH_Y", "NTIMESTEPS"}

def _prompt_missing_fields(system, overrides, skipped_fields, system_names):
    """Interactively ask for any FIELDS_WORTH_ASKING that are still missing.
    Uses InputCollector for numeric fields. Returns (system, overrides, quit_flag)."""
    missing = missing_fields(system, overrides, skipped_fields)
    while missing:
        key, question = missing[0]
        if key == "system" and system_names:
            print(f"\n  Available systems: {', '.join(system_names)}")

        if key in NUMERIC_FIELDS:
            val = InputCollector.get_int(
                f"  {question} ", min_val=1, max_val=100000000, default=None,
            )
            if val is None:
                return system, overrides, True
            overrides[key] = int(val)
        else:
            text, quit = _get_user_input(f"  {question} ")
            if quit:
                return system, overrides, True
            if not text:
                missing.pop(0)
                continue
            if key == "system":
                match = next((s for s in system_names if s.lower() == text.lower()), None)
                if not match:
                    match, extra, _extra_mentioned = parse_prompt_llm(text, system_names)
                    overrides = merge_overrides(overrides, extra)
                system = match or text
            else:
                parts = [p.strip() for p in text.split(",")]
                overrides[key] = [_coerce(p) for p in parts] if len(parts) > 1 else _coerce(text)

        missing = missing_fields(system, overrides, skipped_fields)
    return system, overrides, False


def translate(initial_text, base_dir, interactive=True):
    systems = discover_systems(base_dir)
    system_names = list(systems.keys())

    # --- One-shot mode: no interactive loop ---
    if not interactive:
        regex_system, regex_overrides = parse_prompt_regex(initial_text, system_names, systems=systems)
        llm_system, llm_overrides, _ = parse_prompt_llm(initial_text, system_names)
        system = regex_system or llm_system
        overrides = merge_overrides(regex_overrides, llm_overrides)
        prompt_components = extract_components_from_prompt(initial_text)

        if system and system in systems:
            return {"system": system, "overrides": overrides, "raw_prompt": initial_text}
        if prompt_components:
            overrides["COMPONENTS"] = prompt_components
            overrides["NUMCOMPONENTS"] = len(prompt_components)
            return {"system": "NEW", "overrides": overrides, "raw_prompt": initial_text}
        return {"system": None, "overrides": overrides, "raw_prompt": initial_text}

    # --- Interactive mode: loop until system identified or user quits ---
    prompt = initial_text
    empty_count = 0
    garbage_count = 0

    while True:
        # --- Fast garbage filter: reject obviously invalid input immediately ---
        if classify_input(prompt) == "garbage":
            garbage_count += 1
            if garbage_count >= MAX_GARBAGE_ATTEMPTS:
                print(f"\n  {MAX_GARBAGE_ATTEMPTS} unrecognized inputs in a row. Exiting.")
                return {"system": None, "overrides": {}, "raw_prompt": initial_text}
            _, ai_msg = check_input_quality(prompt)
            remaining = MAX_GARBAGE_ATTEMPTS - garbage_count
            hint = f" ({remaining} attempt{'s' if remaining != 1 else ''} left)" if remaining else ""
            if ai_msg:
                print(f"  {ai_msg}{hint}")
            else:
                print(f"  That doesn't seem related to alloy simulation.{hint}")
                print("  Try something like 'aluminum cobalt' or 'NiAl with 6 phases'.")
            text, quit = _get_user_input()
            if quit:
                return {"system": None, "overrides": {}, "raw_prompt": initial_text}
            prompt = text
            continue

        # --- Input passed the filter: reset garbage counter and proceed ---
        garbage_count = 0

        system, overrides = _try_parse(prompt, system_names, systems)

        # 1. Exact system match found -- clarify missing fields and return
        if system and system in systems:
            system, overrides, quit = _prompt_missing_fields(
                system, overrides, detect_unchanged_fields(prompt), system_names
            )
            if quit:
                return {"system": None, "overrides": {}, "raw_prompt": initial_text}
            return {"system": system, "overrides": overrides, "raw_prompt": initial_text}

        # 2. Elements mentioned but don't match any existing system
        components = extract_components_from_prompt(prompt)
        if components:
            comp_set = set(c.upper() for c in components)
            matched = any(
                comp_set == set(c.upper() for c in systems[n].get("components", []))
                for n in system_names
            )
            if not matched:
                result = _offer_new_system(system_names)
                if result:
                    result["raw_prompt"] = initial_text
                    return result
                # User quit or said no — stop immediately
                return {"system": None, "overrides": {}, "raw_prompt": initial_text}

        # 3. Nothing recognizable -- ask for clarification
        print(f"\n  Available systems: {', '.join(system_names)}")
        text, quit = _get_user_input(
            "  I couldn't match that to an existing system. "
            "Can you describe it differently, or mention a system name? "
        )
        if quit:
            return {"system": None, "overrides": {}, "raw_prompt": initial_text}
        if not text:
            empty_count += 1
            if empty_count >= 3:
                print("  No valid input received. Exiting.")
                return {"system": None, "overrides": {}, "raw_prompt": initial_text}
            continue
        empty_count = 0
        prompt = text


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("prompt", nargs="*", help="One-shot prompt text (omit for interactive mode)")
    ap.add_argument("-o", "--output", default="spec.json", help="Where to write the resulting spec JSON")
    ap.add_argument("--base-dir", default=os.environ.get("MICROSIM_DIR", "/app/microsim_gp"))
    ap.add_argument("--non-interactive", action="store_true",
                     help="Never prompt the user; fill only what parsing finds")
    args = ap.parse_args()

    if not _ollama_available():
        print(f"  [translator] Warning: Ollama not reachable at {OLLAMA_HOST}. "
              f"Falling back to regex-only parsing.")

    if args.prompt:
        text = " ".join(args.prompt)
    else:
        systems = discover_systems(args.base_dir)
        if systems:
            show_systems(systems)
        print("  Describe what you want to simulate, e.g.:")
        print("    'AlZn_eutectic with 6 phases, mesh_x 200, mesh_y 200, timesteps 150000'")
        print("  Or describe a new system: 'aluminum cobalt alloy'")
        print("  Type 'quit' at any time to exit.\n")
        try:
            text = input(">>> ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nGoodbye.")
            sys.exit(0)
        if text.lower() in QUIT_COMMANDS:
            print("Goodbye.")
            sys.exit(0)

    spec = translate(text, args.base_dir, interactive=not args.non_interactive)

    if not spec.get("system"):
        print("\n  No system selected. Exiting.")
        sys.exit(1)

    print("\n  --- Translated spec ---")
    print(f"  System: {spec['system']}")
    if spec["overrides"]:
        for k, v in sorted(spec["overrides"].items()):
            print(f"    {k} = {v}")
    else:
        print("    (no parameter overrides)")

    with open(args.output, "w") as f:
        json.dump(spec, f, indent=2)
    print(f"\n  Wrote spec to {args.output}")


if __name__ == "__main__":
    main()
