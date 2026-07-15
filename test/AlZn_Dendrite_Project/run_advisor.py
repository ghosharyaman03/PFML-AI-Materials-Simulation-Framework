import autogen
import os

# 1. Point AutoGen to your local Ollama instance
config_list = [
    {
        "model": "llama3",                  # Change this to your exact Ollama model name if different (e.g., mistral, phi3)
        "api_base": "http://localhost:11434/v1",  # Local Ollama endpoint port
        "api_type": "openai",
        "api_key": "ollama",                # Ollama doesn't require a real key
    }
]

llm_config = {"config_list": config_list, "cache_seed": None}

# 2. Define the Simulation Orchestrator Agent
orchestrator = autogen.ConversableAgent(
    name="Orchestrator",
    system_message="You are the Simulation Coordinator. You take the user's parameter requests, format them clearly, and pass them to the Materials Expert for validation. Keep responses conversational and brief.",
    llm_config=llm_config,
    human_input_mode="NEVER"
)

# 3. Define the Materials Science Expert Agent
materials_expert = autogen.ConversableAgent(
    name="Materials_Science_Expert",
    system_message="You are an expert in phase-field dendritic growth simulations for Al-Zn alloys. "
                   "When parameters change (like anisotropy strength epsilon or solute distribution), "
                   "you analyze how it impacts tip velocity, segregation, and microstructural arm spacing. Keep responses concise.",
    llm_config=llm_config,
    human_input_mode="NEVER"
)

# 4. Get the modification parameter from the user dynamically
print("\n" + "="*50)
print("   AutoGen + Ollama Al-Zn Advisor Loaded")
print("="*50)
user_tweak = input("\nEnter your Al-Zn parameter changes (e.g., Increase anisotropy strength to 0.05): ")

simulation_context = f"""
Workspace Context: Modifying phase-field microstructural evolution parameters.
Alloy System: Al-Zn (Aluminum-Zinc) dendritic solidification.
Proposed Change: {user_tweak}
"""

print("\n--- Starting Local AutoGen Orchestrator Session ---\n")

# 5. Kick off the agent-to-agent conversation
orchestrator.initiate_chat(
    recipient=materials_expert,
    message=f"Please analyze this simulation tweak based on our files:\n{simulation_context}",
    max_turns=2
)