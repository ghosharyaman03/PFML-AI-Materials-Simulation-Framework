
def _try_parse(prompt, system_names, systems):
    """Parse a prompt through regex + LLM. Returns (system, overrides,
    mentioned_no_value) -- the third item lists parameters the LLM
    recognized as referenced (via its baked MicroSim domain knowledge)
    but couldn't ground a specific number for, e.g. "interface relaxes
    faster" -> tau with no number attached."""
    regex_system, regex_overrides = parse_prompt_regex(prompt, system_names, systems=systems)
    llm_system, llm_overrides, mentioned_no_value = parse_prompt_llm(prompt, system_names)

    # General safety net (not tied to any one system): regex's own
    # matching already enforces component consistency, including the
    # SYSTEM_NAME_ALIASES exception for known-mislabeled presets like
    # NiAlCr -- but the LLM has no such check at all, so it can
    # confidently hallucinate a real system name that doesn't actually
    # contain what the user mentioned (e.g. "copper zinc" -> AlZn_dendrite,
    # which has no copper). Only applies when regex didn't already find a
    # system -- if it did, that path is already trustworthy.
    if llm_system and not regex_system and systems and llm_system in systems:
        mentioned = {
            sym.upper() for word, sym in ELEMENT_KEYWORDS.items()
            if re.search(rf"\b{word}\b", prompt, re.IGNORECASE)
        }
        # Also catch short symbols typed directly (e.g. "ti nb"), not just
        # full element names -- otherwise this safety net can't see them
        # and an LLM hallucination like "ti nb" -> NiNb slips through.
        mentioned |= {
            sym for sym in extract_components_from_prompt(prompt)
        }
        mentioned = {m.upper() for m in mentioned}
        if mentioned:
            comps = systems[llm_system]["config"].get("COMPONENTS", [])
            if isinstance(comps, str):
                comps = [comps]
            comp_set = {str(c).upper() for c in comps}
            if comp_set and not mentioned.issubset(comp_set):
                print(f"  [translator] Discarding AI-suggested system '{llm_system}' -- "
                      f"you mentioned {sorted(mentioned)} but that system's actual "
                      f"components are {sorted(comp_set)}.")
                llm_system = None

    system = regex_system or llm_system
    overrides = merge_overrides(regex_overrides, llm_overrides)
    mentioned_no_value = [k for k in mentioned_no_value if k not in overrides]
    return system, overrides, mentioned_no_value


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

        system, overrides, mentioned_no_value = _try_parse(prompt, system_names, systems)

        # 1. Exact system match found -- ask about anything the model
        #    recognized as referenced but couldn't ground a number for
        #    (e.g. "interface relaxes faster" -> tau), then clarify the
        #    always-required fields, then return.
        if system and system in systems:
            for key in mentioned_no_value:
                if key in overrides:
                    continue  # resolved by regex or a later grounded value
                text, quit = _get_user_input(
                    f"  You mentioned something related to '{key}' -- "
                    f"what value would you like? (comma-separate multiple "
                    f"values, or leave blank to skip): "
                )
                if quit:
                    return {"system": None, "overrides": {}, "raw_prompt": initial_text}
                if not text:
                    continue
                parts = [p.strip() for p in text.split(",") if p.strip()]
                overrides[key] = [_coerce(p) for p in parts] if len(parts) > 1 else _coerce(parts[0])

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
