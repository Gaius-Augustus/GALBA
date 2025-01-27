import os
import sys
import shutil
from pathlib import Path
import logging

def locate_augustus_config(working_dir):
    """
    Attempt to find and/or create a writable AUGUSTUS_CONFIG_PATH.
    Returns the final config path that is guaranteed to be writable,
    and sets os.environ["AUGUSTUS_CONFIG_PATH"] to that path.

    Logic summary:
      1) If $AUGUSTUS_CONFIG_PATH is set and is a directory, try that first.
      2) If not set or not valid, guess from an installed 'augustus' or
         from a known default location.
      3) Check if the config path is writable.
         - If not writable, copy it to e.g. $HOME/.augustus or
           working_dir/.augustus, then point $AUGUSTUS_CONFIG_PATH to the copy.
    """
    logger = logging.getLogger("galba_pipeline.augustus")

    # 1) Check environment variable
    aug_path = os.environ.get("AUGUSTUS_CONFIG_PATH", "")
    if aug_path:
        logger.info(f"Found $AUGUSTUS_CONFIG_PATH in environment: {aug_path}")
    else:
        # 2) "Guess" a location if environment is missing
        #    e.g. if 'augustus' is installed at /usr/local/bin/augustus,
        #    the config might be at /usr/local/bin/../config or /usr/share/augustus/config
        # For brevity, we’ll pick a common fallback:
        guess_path = "/usr/share/augustus/config"
        logger.info(
            f"Environment variable $AUGUSTUS_CONFIG_PATH not set. "
            f"Guessing default: {guess_path}"
        )
        aug_path = guess_path

    # Ensure aug_path is a real absolute directory
    aug_path = os.path.abspath(aug_path)
    if not os.path.isdir(aug_path):
        logger.warning(
            f"Guessed or environment config path '{aug_path}' is not a directory. "
            "We may need to do further checks or fallback to other guesses."
        )
        # Possibly fallback to a second guess, or fail
        # For now, let’s just create an empty directory to avoid immediate crash:
        # (In practice, you might try a second guess or exit with an error.)
        # return or raise an error:
        logger.error(f"No valid AUGUSTUS config folder found at {aug_path}.")
        sys.exit(1)

    # 3) Check if it's writable:
    if is_writable(aug_path):
        logger.info(f"AUGUSTUS_CONFIG_PATH is writable: {aug_path}")
        # We can keep this location
        os.environ["AUGUSTUS_CONFIG_PATH"] = aug_path
        return aug_path
    else:
        logger.warning(f"AUGUSTUS config path {aug_path} is not writable.")
        # 4) Attempt to copy it to a writable location (e.g. home or local)
        # Option A: user’s home directory
        home_dir = Path.home()
        alt_path = home_dir / ".augustus"
        # If that fails for some reason (read-only container?), you can try `working_dir/.augustus`
        # or anything else that is definitely writable, e.g.:
        # alt_path = Path(working_dir) / ".augustus"

        if alt_path.exists():
            # If .augustus already exists, we might check if it’s the same config or we forcibly remove it
            logger.info(f"Found existing {alt_path}. We will not overwrite it unless it’s stale.")
        else:
            logger.info(f"Copying {aug_path} => {alt_path} to get a writable config folder.")
            try:
                shutil.copytree(aug_path, alt_path)
            except Exception as e:
                logger.error(
                    f"Failed to copy {aug_path} to {alt_path}.\n"
                    "Try specifying a different location or set $AUGUSTUS_CONFIG_PATH to something writable."
                )
                sys.exit(1)

        # Now we check that alt_path is definitely writable
        if not is_writable(str(alt_path)):
            logger.error(
                f"Copied config to {alt_path}, but that directory is still not writable!"
            )
            sys.exit(1)

        # 5) If success, update environment
        os.environ["AUGUSTUS_CONFIG_PATH"] = str(alt_path)
        logger.info(
            f"Reset $AUGUSTUS_CONFIG_PATH to {alt_path} because the original was not writable."
        )
        return str(alt_path)

def is_writable(path):
    """
    Check if 'path' is a directory and if the current user has write + execute
    permissions on it (common requirement to create files/subdirectories).
    Returns True if yes, otherwise False.
    Fragile, there are some file systems where this may not work!
    """
    if not os.path.isdir(path):
        return False
    
    # 'write' permission alone is not enough for creating files within a directory;
    # you typically also need 'execute' permission on the directory.
    can_write = os.access(path, os.W_OK)
    can_execute = os.access(path, os.X_OK)

    return can_write and can_execute