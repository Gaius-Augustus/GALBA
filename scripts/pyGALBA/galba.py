#!/usr/bin/env python3

"""
main.py - Entry point for the refactored GALBA pipeline.
"""

import logging
import sys
import os

from pipeline.arguments import parse_arguments
from pipeline.dependencies import check_dependencies
from pipeline.augustus_config import locate_augustus_config
from pipeline.fileformats import check_fileformats
from pipeline.checkpoint import load_state, save_state
from pipeline.runminiprot import run_miniprot
from pipeline.hints import generate_hints
from pipeline.training import run_training
from pipeline.prediction import run_prediction
from pipeline.postprocessing import postprocess_predictions

def main():
    # 1) Parse command-line arguments
    args = parse_arguments(sys.argv[1:])
    check_dependencies(args) # checks for required executables and scripts
    check_fileformats(args) # checks for required file formats
    # 2) Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s'
    )
    logger = logging.getLogger("galba_pipeline")

    # 3) Create working directory
    try:
        if not os.path.exists(args.workingdir):
            os.makedirs(args.workingdir)
    except Exception as e:
        logger.error(f"Error creating working directory: {e}")
        sys.exit(1)
    # 3.1) Locate AUGUSTUS config (possibly copying into workingdir)
    locate_augustus_config(args.workingdir)

    # 4) Load or create pipeline state to allow resuming
    pipeline_state = load_state(args.workingdir)
    
    # 5) Potentially skip steps if checkpoint says they're done
    #    e.g., if pipeline_state["alignment_done"] is True, skip alignment, etc.
    
    # Step A: Protein alignment (only if not done)
    if not pipeline_state.get("runminiprot_done", False):
        logger.info("Starting runminiprot step...")
        run_miniprot(args)
        pipeline_state["runminiprot_done"] = True
        save_state(args.workingdir, pipeline_state)
    
    # Step B: Generate/merge hints
    if not pipeline_state.get("hints_done", False):
        logger.info("Generating hints...")
        generate_hints(args)
        pipeline_state["hints_done"] = True
        save_state(args.workingdir, pipeline_state)
    
    # Step C: Training
    if not pipeline_state.get("training_done", False):
        logger.info("Running training...")
        run_training(args)
        pipeline_state["training_done"] = True
        save_state(args.workingdir, pipeline_state)
    
    # Step D: Prediction
    if not pipeline_state.get("prediction_done", False):
        logger.info("Running final prediction...")
        run_prediction(args)
        pipeline_state["prediction_done"] = True
        save_state(args.workingdir, pipeline_state)
    
    # Step E: Post-processing
    if not pipeline_state.get("postprocessing_done", False):
        logger.info("Post-processing predictions...")
        postprocess_predictions(args)
        pipeline_state["postprocessing_done"] = True
        save_state(args.workingdir, pipeline_state)

    logger.info("GALBA pipeline finished successfully!")

if __name__ == "__main__":
    main()
