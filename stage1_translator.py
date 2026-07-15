#!/usr/bin/env python3
"""
stage1_translator.py
=====================
Turns a user's free-form paragraph into a structured simulation spec, e.g.

  "hey can u run something with that aluminum zinc alloy, the eutectic one,
   give it like 6 phases and make the grid a bit bigger, maybe 200 by 200,
   and let it run for a good while, 150000 steps"

      -->

  {"system": "AlZn_eutectic",
   "overrides": {"NUMPHASES": 6, "MESH_X": 200, "MESH_Y": 200, "NTIMESTEPS": 150000}}

It never invents thermodynamic values, TDB choices, or anything belonging
to stage2/stage3 -- it only figures out *which* system and *which*
parameters the user wants changed. If it can't confidently fill in the
handful of fields that matter (system, NUMPHASES, MESH_X, MESH_Y,
NTIMESTEPS), it asks the user directly instead of guessing.

Uses a local Ollama model (small, e.g. 2B) to help parse messy language,
with the existing regex parser (shared_core.parse_prompt-equivalent) as
a free, fast first pass and a safety net if Ollama is unavailable.

Usage:
  Interactive:  python3 stage1_translator.py
  One-shot:     python3 stage1_translator.py "AlZn_eutectic with 6 phases..." -o spec.json

Env vars:
  MICROSIM_DIR    Path to microsim_gp directory (default: /app/microsim_gp)
  OLLAMA_HOST     Base URL of the Ollama server (default: http://localhost:11434)
  OLLAMA_MODEL    Model to use for parsing (default: gemma2:2b)
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

# The handful of fields worth actively asking about if we can't infer
# them. Everything else in PARAM_ALIASES/SECTIONS can still be picked up
# by regex/LLM parsing, but we don't interrogate the user over minor
# knobs -- only the ones that meaningfully change what gets simulated.
FIELDS_WORTH_ASKING = [
    ("system", "Which alloy system do you want to simulate?"),
    ("NUMPHASES", "How many phases should the simulation have?"),
    ("MESH_X", "What grid/mesh width (MESH_X) do you want?"),
    ("MESH_Y", "What grid/mesh height (MESH_Y) do you want?"),
    ("NTIMESTEPS", "How many timesteps should it run for?"),
]

# Element-name keywords, used as a fallback when no preset name literally
# appears in the prompt (e.g. "aluminum copper" instead of "Model_Solidification").
ELEMENT_KEYWORDS = {
    "aluminum": "Al", "aluminium": "Al",
    "zinc": "Zn",
    "nickel": "Ni",
    "copper": "Cu",
    "molybdenum": "Mo",
    "niobium": "Nb",
    "chromium": "Cr",
}
SYSTEM_NAME_ALIASES = {
    "NiAlCr": ["chromium", "chromate", "cr alloy"],
}
NUMBER_WORDS = {
    "one": 1, "two": 2, "three": 3, "four": 4, "five": 5,
    "six": 6, "seven": 7, "eight": 8, "nine": 9, "ten": 10,
}
# Generic naming tokens that shouldn't be used to disambiguate between
# presets sharing the same components (e.g. "Model" in "Model_Solidification"
# vs "Model_Precipitation" -- too common to mean anything on its own).
GENERIC_NAME_TOKENS = {"model", "example", "examples", "system", "alloy"}

# Phrases indicating the user wants a field left at the preset default,
# e.g. "keep phases the same" -- treated as answered, not missing.
UNCHANGED_SIGNAL_WORDS = {
    "NUMPHASES": ["phase", "phases"],
    "MESH_X": ["mesh"],
    "MESH_Y": ["mesh"],
    "NTIMESTEPS": ["timestep", "timesteps", "step", "steps"],
}
UNCHANGED_PATTERN = r"(same|unchanged|as[- ]is|no change|keep it)"


def detect_unchanged_fields(prompt):
    """Return the set of FIELDS_WORTH_ASKING keys the user explicitly
    said to leave at the preset default (e.g. "keep phases the same"),
    so translate() doesn't ask about them even though no number was given."""
    skipped = set()
    for key, topic_words in UNCHANGED_SIGNAL_WORDS.items():
        for word in topic_words:
            if re.search(rf"\b{word}\b[^.]{{0,25}}\b{UNCHANGED_PATTERN}\b", prompt, re.IGNORECASE) or \
               re.search(rf"\b{UNCHANGED_PATTERN}\b[^.]{{0,25}}\b{word}\b", prompt, re.IGNORECASE):
                skipped.add(key)
                break
    return skipped
# ============================================================================
# Regex fallback parser (fast, free, no dependencies)
# ============================================================================

def parse_prompt_regex(prompt, system_names=None, systems=None):
    """Extract a system name and numeric overrides from natural language
    using plain regex/alias matching. Deliberately conservative -- this
    is the safety net under the LLM pass, not the primary parser.
    Returns (system_name | None, overrides_dict).
    """
    system_names = system_names or []
    preset = None
    overrides = {}

    for name in sorted(system_names, key=len, reverse=True):
        pattern = re.escape(name).replace("_", r"[_ ]?")
        if re.search(rf"\b{pattern}\b", prompt, re.IGNORECASE):
            preset = name
            break

    # Fallback: no preset name was literally mentioned -- try inferring
    # it from element keywords (e.g. "aluminum copper" -> Al, Cu), then
    # disambiguate between same-component presets using a distinctive
    # word from each candidate's own name (e.g. "solidification").
    if not preset:
        for name, aliases in SYSTEM_NAME_ALIASES.items():
            if name in system_names and any(
                re.search(rf"\b{re.escape(alias)}\b", prompt, re.IGNORECASE) for alias in aliases
            ):
                preset = name
                break

    # Manual name aliases -- for presets whose stored COMPONENTS don't
    # match what their name implies, so element-keyword matching alone
    # would never resolve them.
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
    """Single-shot call to a local Ollama model. system_prompt is optional
    -- the production model (pfml-parser) already has the valid systems/
    params baked into its own Modelfile SYSTEM block via
    build_modelfile.py, so omitting it here lets Ollama use that baked
    knowledge. Passing system_prompt explicitly overrides it (useful for
    testing against a raw base model that wasn't built with build_modelfile.py)."""
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


def parse_prompt_llm(prompt, system_names, known_overrides=None):
    """Ask the local Ollama model to extract system + overrides.
    Returns (system_name | None, overrides_dict) or (None, {}) on failure.
    Never invents values for fields the user didn't mention.
    """
    if not _ollama_available():
        return None, {}

    all_params = set()
    for _, keys in SECTIONS:
        all_params.update(keys)
    all_params.update(PARAM_ALIASES.values())

    raw = ask_ollama(prompt)
    if not raw:
        return None, {}

    try:
        text = re.sub(r"```json\s*|```\s*", "", raw).strip()
        parsed = json.loads(text)
    except Exception:
        return None, {}

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

    return resolved_system, coerced


# ============================================================================
# Interactive clarification loop
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


def translate(initial_text, base_dir, interactive=True):
    """Run the translate-and-clarify loop. Returns a spec dict:
    {"system": ..., "overrides": {...}, "raw_prompt": ...}
    """
    systems = discover_systems(base_dir)
    system_names = list(systems.keys())

    regex_system, regex_overrides = parse_prompt_regex(initial_text, system_names, systems=systems)
    llm_system, llm_overrides = parse_prompt_llm(initial_text, system_names)

    system = regex_system or llm_system
    overrides = merge_overrides(regex_overrides, llm_overrides)
    skipped_fields = detect_unchanged_fields(initial_text)

    if not interactive:
        return {"system": system, "overrides": overrides, "raw_prompt": initial_text}

    # --- Clarify anything important that's still missing ---
    missing = missing_fields(system, overrides, skipped_fields)
    while missing:
        key, question = missing[0]
        if key == "system" and system_names:
            print(f"\n  Available systems: {', '.join(system_names)}")
        try:
            answer = input(f"  {question} ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\n  Stopping -- using what's been gathered so far.")
            break

        if not answer:
            # User declined to specify -- stop asking about this field
            # and fall back to whatever the preset default provides.
            missing.pop(0)
            continue

        if key == "system":
            match = next((s for s in system_names if s.lower() == answer.lower()), None)
            if not match:
                # Try the LLM once more with the direct answer
                match, extra = parse_prompt_llm(answer, system_names)
                overrides = merge_overrides(overrides, extra)
            system = match or answer
        else:
            parts = [p.strip() for p in answer.split(",")]
            overrides[key] = [_coerce(p) for p in parts] if len(parts) > 1 else _coerce(answer)

        missing = missing_fields(system, overrides, skipped_fields)

    return {"system": system, "overrides": overrides, "raw_prompt": initial_text}


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
        try:
            text = input(">>> ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nGoodbye.")
            sys.exit(0)

    spec = translate(text, args.base_dir, interactive=not args.non_interactive)

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
