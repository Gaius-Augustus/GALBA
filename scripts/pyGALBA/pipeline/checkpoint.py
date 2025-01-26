# pipeline/checkpoint.py
import os
import json

def load_state(workingdir):
    """
    Load pipeline state (JSON) from disk if it exists; else return empty dict.
    """
    state_file = os.path.join(workingdir, "pipeline_state.json")
    if os.path.exists(state_file):
        with open(state_file, "r") as f:
            return json.load(f)
    return {}

def save_state(workingdir, state_dict):
    """
    Write the pipeline state (JSON) to disk.
    """
    state_file = os.path.join(workingdir, "pipeline_state.json")
    with open(state_file, "w") as f:
        json.dump(state_dict, f, indent=2)
