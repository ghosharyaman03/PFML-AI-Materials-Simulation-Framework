import autogen
import re

# =====================================================================
# STRONGHOLD 1: EXPLICIT IN-MEMORY METALLURGICAL REGISTRY (RAW NUMBERS)
# =====================================================================
WORKSPACE_REGISTRY = {
    "alzn": {
        "system_name": "Al-Zn (Aluminum-Zinc) Dendrite Simulation Workspace",
        "mesh": {"dimensions": [200, 100, 250], "delta_x": 2e-8, "delta_y": 2e-8, "delta_z": 2e-8},
        "discretization": {"delta_t": 2e-8, "ntimesteps": 50000, "savet": 10000},
        "kinetics": {"filling_temperature": 857, "equilibrium_temperature": 870, "under_cooling": 13},
        "equilibrium_compositions": {
            "alpha_solidus": 0.926,
            "liquid_liquidus": 0.817
        },
        "thermodynamics": {"slopes_solidus_alpha": 772.34, "slopes_liquidus_alpha": 299.26, "diffusivity_liquid": 1e-9},
        "potential_parameters": {"A_alpha": 30015.44, "A_liquid": 12808.21, "function_F": 3, "function_W": 1},
        "anisotropy": {"function_anisotropy": 1, "anisotropy_type": 4, "dab": 0.02, "rotation_angle_degrees": 45}
    },
    "nial": {
        "system_name": "Ni-Al (Nickel-Aluminum) Isothermal Workspace",
        "mesh": {"dimensions": [100, 100, 1], "delta_x": 1e-7, "delta_y": 1e-7, "delta_z": 1e-7},
        "discretization": {"delta_t": 1e-7, "ntimesteps": 20000, "savet": 1000},
        "kinetics": {"filling_temperature": 1717, "equilibrium_temperature": 1718, "under_cooling": 1},
        "geometry": "FILLCYLINDER = {0,0,0,0,0,5} seed at origin",
        "thermodynamics": {
            "1717.0": {"Al_fcc": 0.0985458, "Al_liq": 0.107529, "hsn_fcc": 315006, "hsn_liq": 275167},
            "1717.1": {"Al_fcc": 0.0981515, "Al_liq": 0.107076, "hsn_fcc": 314957, "hsn_liq": 275160},
            "1717.2": {"Al_fcc": 0.0977544, "Al_liq": 0.106621, "hsn_fcc": 314912, "hsn_liq": 275157},
            "1717.3": {"Al_fcc": 0.0973546, "Al_liq": 0.106163, "hsn_fcc": 314871, "hsn_liq": 275159},
            "1717.4": {"Al_fcc": 0.0969520, "Al_liq": 0.105702, "hsn_fcc": 314834, "hsn_liq": 275164},
            "1717.5": {"Al_fcc": 0.0965467, "Al_liq": 0.105237, "hsn_fcc": 314800, "hsn_liq": 275175},
            "1717.6": {"Al_fcc": 0.0961384, "Al_liq": 0.104769, "hsn_fcc": 314771, "hsn_liq": 275190},
            "1717.7": {"Al_fcc": 0.0957272, "Al_liq": 0.104299, "hsn_fcc": 314747, "hsn_liq": 275210},
            "1717.8": {"Al_fcc": 0.0953131, "Al_liq": 0.103824, "hsn_fcc": 314727, "hsn_liq": 275235},
            "1717.9": {"Al_fcc": 0.0948960, "Al_liq": 0.103347, "hsn_fcc": 314712, "hsn_liq": 275265},
            "1718.0": {"Al_fcc": 0.0944758, "Al_liq": 0.102866, "hsn_fcc": 314701, "hsn_liq": 275301}
        }
    },
    "ninb": {
        "system_name": "Ni-Nb (Nickel-Niobium) High-Fidelity Simulation Workspace",
        "mesh": {"dimensions": [200, 200, 1], "delta_x": 3e-8},
        "discretization": {"delta_t": 162e-9, "ntimesteps": 2000000},
        "thermodynamics": {
            "1593": {"Nb_fcc": 0.0996303, "Nb_liq": 0.123211, "hsn_fcc": 596363, "hsn_liq": 549773},
            "1594": {"Nb_fcc": 0.0989523, "Nb_liq": 0.122454, "hsn_fcc": 597743, "hsn_liq": 550782},
            "1595": {"Nb_fcc": 0.0982729, "Nb_liq": 0.121695, "hsn_fcc": 599139, "hsn_liq": 551803},
            "1596": {"Nb_fcc": 0.0975920, "Nb_liq": 0.120935, "hsn_fcc": 600552, "hsn_liq": 552835},
            "1597": {"Nb_fcc": 0.0969096, "Nb_liq": 0.120173, "hsn_fcc": 601982, "hsn_liq": 553880},
            "1598": {"Nb_fcc": 0.0962259, "Nb_liq": 0.119409, "hsn_fcc": 603429, "hsn_liq": 554938},
            "1599": {"Nb_fcc": 0.0955406, "Nb_liq": 0.118643, "hsn_fcc": 604894, "hsn_liq": 556009},
            "1600": {"Nb_fcc": 0.0948540, "Nb_liq": 0.117876, "hsn_fcc": 606378, "hsn_liq": 557092},
            "1601": {"Nb_fcc": 0.0941658, "Nb_liq": 0.117107, "hsn_fcc": 607880, "hsn_liq": 558189},
            "1602": {"Nb_fcc": 0.0934763, "Nb_liq": 0.116336, "hsn_fcc": 609402, "hsn_liq": 559300},
            "1603": {"Nb_fcc": 0.0927853, "Nb_liq": 0.115563, "hsn_fcc": 610943, "hsn_liq": 560425},
            "1604": {"Nb_fcc": 0.0920928, "Nb_liq": 0.114788, "hsn_fcc": 612504, "hsn_liq": 561565},
            "1605": {"Nb_fcc": 0.0913989, "Nb_liq": 0.114012, "hsn_fcc": 614087, "hsn_liq": 562719},
            "1606": {"Nb_fcc": 0.0907035, "Nb_liq": 0.113234, "hsn_fcc": 615690, "hsn_liq": 563888},
            "1607": {"Nb_fcc": 0.0900067, "Nb_liq": 0.112453, "hsn_fcc": 617315, "hsn_liq": 565072},
            "1608": {"Nb_fcc": 0.0893084, "Nb_liq": 0.111672, "hsn_fcc": 618962, "hsn_liq": 566272},
            "1609": {"Nb_fcc": 0.0886087, "Nb_liq": 0.110888, "hsn_fcc": 620633, "hsn_liq": 567488},
            "1610": {"Nb_fcc": 0.0879075, "Nb_liq": 0.110102, "hsn_fcc": 622326, "hsn_liq": 568720},
            "1611": {"Nb_fcc": 0.0872049, "Nb_liq": 0.109315, "hsn_fcc": 624044, "hsn_liq": 569970},
            "1612": {"Nb_fcc": 0.0865008, "Nb_liq": 0.108526, "hsn_fcc": 625786, "hsn_liq": 571236},
            "1613": {"Nb_fcc": 0.0857953, "Nb_liq": 0.107735, "hsn_fcc": 627554, "hsn_liq": 572520},
            "1614": {"Nb_fcc": 0.0850883, "Nb_liq": 0.106942, "hsn_fcc": 629348, "hsn_liq": 57322},
            "1615": {"Nb_fcc": 0.0843799, "Nb_liq": 0.106147, "fcc": 631168, "hsn_liq": 575143},
            "1616": {"Nb_fcc": 0.0836700, "Nb_liq": 0.105351, "hsn_fcc": 633017, "hsn_liq": 576482},
            "1617": {"Nb_fcc": 0.0829587, "Nb_liq": 0.104552, "hsn_fcc": 634893, "hsn_liq": 577841},
            "1618": {"Nb_fcc": 0.0822459, "Nb_liq": 0.103752, "hsn_fcc": 636798, "hsn_liq": 579220},
            "1619": {"Nb_fcc": 0.0815318, "Nb_liq": 0.102950, "hsn_fcc": 638733, "hsn_liq": 580619},
            "1620": {"Nb_fcc": 0.0808161, "Nb_liq": 0.102145, "hsn_fcc": 640699, "hsn_liq": 582040},
            "1621": {"Nb_fcc": 0.0800991, "Nb_liq": 0.101339, "hsn_fcc": 642696, "hsn_liq": 583481},
            "1622": {"Nb_fcc": 0.0793806, "Nb_liq": 0.100532, "hsn_fcc": 644726, "hsn_liq": 584945}
        }
    }
}

# =====================================================================
# STRONGHOLD 2: EXACT ARCHITECTURAL INJECTION PARAMETERS
# =====================================================================
config_list = [{
    'model': 'gemma2:2b',
    'base_url': 'http://host.docker.internal:11434/v1',
    'api_key': 'ollama',
}]

llm_config = {
    'config_list': config_list,
    'temperature': 0.0, 
    'timeout': 600,
}

user_proxy = autogen.UserProxyAgent(
    name='user_proxy',
    human_input_mode='NEVER',
    max_consecutive_auto_reply=2,  
    is_termination_msg=lambda x: "exitcode: 0" in x.get("content", "") or "FINAL_ANSWER" in x.get("content", ""),
    code_execution_config={
        "work_dir": "/workspace/test",
        "use_docker": False        
    }
)

def enforce_filename_hook(recipient, messages, sender, config):
    if messages and messages[-1].get("content") and hasattr(recipient, 'enforced_filename'):
        target = recipient.enforced_filename
        messages[-1]["content"] = re.sub(
            r"['\"][\w-]+\.(?:in|txt|cfg|conf)['\"]", 
            f"'{target}'", 
            messages[-1]["content"]
        )
    return False, None

user_proxy.register_reply([autogen.Agent, None], reply_func=enforce_filename_hook, position=0)

# =====================================================================
# 3. RUNTIME SYSTEM EXECUTION LOOP
# =====================================================================
print('='*70)
print('[System] Hardcoded Memory Matrix Engaged. Pipeline Active.')
print('='*70)

user_input = input("\nTask / Question Statement: ")

if user_input.strip():
    filename_match = re.search(r'([\w-]+\.(?:in|txt|cfg|conf))', user_input)
    enforced_filename = filename_match.group(1) if filename_match else 'simulation.in'
    
    user_proxy.enforced_filename = enforced_filename

    # Explicit string construction completely eliminates triple-quote escaping bugs
    base_system_message = (
        "You are a dual-purpose simulation analyst and configuration generator.\n\n"
        "TASK 1: MATHEMATICAL DELTAS\n"
        "If the user asks to calculate or compute a change:\n"
        "1. Print the calculation steps and append 'FINAL_ANSWER' at the end. Stop.\n\n"
        "TASK 2: CONFIGURATION FILE CREATION\n"
        "1. Locate the requested alloy parameters inside the database registry context.\n"
        "2. You MUST write out EVERY SINGLE structural parameter group belonging to that alloy.\n"
        "3. Provide a clean Python code block that explicitly maps all these values.\n"
        "4. Use exactly '__TARGET_FILE__' as your file output name.\n\n"
        "Your response MUST contain an executable Python segment following this layout precisely:\n\n"
        "```python\n"
        "config_data = \"\"\"[SYSTEM]\n"
        "system_name = ...\n\n"
        "[MESH]\n"
        "dimensions = ...\n"
        "delta_x = ...\n"
        "delta_y = ...\n"
        "delta_z = ...\n\n"
        "[DISCRETIZATION]\n"
        "delta_t = ...\n"
        "ntimesteps = ...\n"
        "savet = ...\n\n"
        "[KINETICS]\n"
        "filling_temperature = ...\n"
        "equilibrium_temperature = ...\n"
        "under_cooling = ...\n\n"
        "[GEOMETRY]\n"
        "layout = ...\n"
        "\"\"\"\n"
        "with open('__TARGET_FILE__', 'w') as f:\n"
        "    f.write(config_data)\n"
        "print(\"SUCCESS: Input file generated.\")\n"
        "```\n"
        "5. Do not provide conversational filler outside of your code block."
    )

    orchestrator = autogen.AssistantAgent(
        name='orchestrator',
        llm_config=llm_config,
        system_message=base_system_message.replace('__TARGET_FILE__', enforced_filename)
    )

    universal_runtime_prompt = f"""
    CONTEXT DATA MATRIX REGISTRY:
    {WORKSPACE_REGISTRY}
    
    USER COMPUTE REQUEST:
    {user_input}
    """
    
    print(f'\n[System] Parsing data grid layers... Target Filename Enforced: {enforced_filename}\n')
    user_proxy.initiate_chat(orchestrator, message=universal_runtime_prompt)
else:
    print("[System] Input error. Workspace query is empty.")