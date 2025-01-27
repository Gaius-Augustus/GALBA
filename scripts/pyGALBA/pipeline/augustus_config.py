import os
import sys
import shutil
from pathlib import Path
import logging
import urllib.request

def locate_augustus_config(working_dir):
    """
    Attempt to find and/or create a writable AUGUSTUS_CONFIG_PATH.
    Returns the final config path that is guaranteed to be writable,
    and sets os.environ["AUGUSTUS_CONFIG_PATH"] to that path.

    Logic summary:
      1) If $AUGUSTUS_CONFIG_PATH is set and is a directory, use it first.
      2) If not set or invalid, guess from e.g. /usr/share/augustus/config.
      3) Check if the config path is writable:
         - If not writable, attempt to copy it into:
             a) ~/.augustus
             b) <working_dir>/.augustus
         in that order, until success. Then point $AUGUSTUS_CONFIG_PATH to that copy.
      4) Finally, ensure that aug_cmdln_parameters.json is present in the config,
         and if missing, attempt to download it from GitHub.
    """

    logger = logging.getLogger("galba_pipeline.augustus")

    # 1) Check environment variable
    aug_path = os.environ.get("AUGUSTUS_CONFIG_PATH", "")
    if aug_path:
        logger.info(f"Found $AUGUSTUS_CONFIG_PATH in environment: {aug_path}")
    else:
        # 2) "Guess" a location if environment is missing
        guess_path = "/usr/share/augustus/config"
        logger.info(
            f"Environment variable $AUGUSTUS_CONFIG_PATH not set. "
            f"Guessing default: {guess_path}"
        )
        aug_path = guess_path

    aug_path = os.path.abspath(aug_path)
    if not os.path.isdir(aug_path):
        logger.error(
            f"Specified or guessed AUGUSTUS config path '{aug_path}' "
            "is not a valid directory. Aborting."
        )
        sys.exit(1)

    # 3) Check if the config path is writable:
    if is_writable(aug_path):
        logger.info(f"AUGUSTUS_CONFIG_PATH is writable: {aug_path}")
        os.environ["AUGUSTUS_CONFIG_PATH"] = aug_path
        ensure_aug_cmdln_params_exists(Path(aug_path))
        return aug_path

    # If not writable, define fallback directories in PREFERRED order:
    fallback_dirs = [
        Path.home() / ".augustus",
        Path(working_dir) / ".augustus"
    ]
    success_path = None

    for alt_path in fallback_dirs:
        logger.info(f"Attempting to copy config {aug_path} to {alt_path} for write access.")
        if attempt_copy_and_set(aug_path, alt_path):
            success_path = alt_path
            break

    if not success_path:
        logger.error("Failed to copy AUGUSTUS config to any fallback location!")
        sys.exit(1)

    # success_path is now guaranteed writable
    ensure_aug_cmdln_params_exists(success_path)
    return str(success_path)

def attempt_copy_and_set(original_path, alt_path):
    """
    Copy 'original_path' to 'alt_path' if needed, then check if 'alt_path' is writable.
    If everything works, set $AUGUSTUS_CONFIG_PATH to alt_path.
    Return True on success, False otherwise.
    """
    logger = logging.getLogger("galba_pipeline.augustus")

    if alt_path.exists():
        logger.info(f"{alt_path} already exists; not overwriting it.")
    else:
        logger.info(f"Copying {original_path} => {alt_path}")
        try:
            shutil.copytree(original_path, alt_path)
        except Exception as e:
            logger.warning(
                f"Could not copy {original_path} to {alt_path}.\n"
                f"Reason: {e}"
            )
            return False

    if not is_writable(str(alt_path)):
        logger.warning(f"Copied to {alt_path}, but that directory is not writable.")
        return False

    os.environ["AUGUSTUS_CONFIG_PATH"] = str(alt_path)
    logger.info(
        f"Reset $AUGUSTUS_CONFIG_PATH to {alt_path} because original was not writable."
    )
    return True

def is_writable(path):
    """
    Check if 'path' is a directory and if the current user has write + execute
    permissions on it. Returns True if yes, otherwise False.
    """
    if not os.path.isdir(path):
        return False
    can_write = os.access(path, os.W_OK)
    can_execute = os.access(path, os.X_OK)
    return can_write and can_execute

def ensure_aug_cmdln_params_exists(alt_path):
    """
    If 'aug_cmdln_parameters.json' is missing from alt_path/parameters/,
    attempt to download it from a known URL, place it there, and log an info message.
    """
    logger = logging.getLogger("galba_pipeline.augustus")
    param_file = alt_path / "parameters" / "aug_cmdln_parameters.json"
    
    if param_file.is_file():
        logger.info(f"AUGUSTUS param file is present: {param_file}")
        return

    # If missing, create parent folder (if needed) and download the file
    logger.warning(
        f"AUGUSTUS param file {param_file} is missing. "
        "We will attempt to download it from GitHub now."
    )
    param_file.parent.mkdir(parents=True, exist_ok=True)

    config_url = (
        "https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/config/parameters/aug_cmdln_parameters.json"
    )

    try:
        logger.info(f"Downloading {config_url} to {param_file}...")
        urllib.request.urlretrieve(config_url, param_file)
        logger.info("Successfully downloaded aug_cmdln_parameters.json.")
    except Exception as e:
        logger.error(
            f"Failed to download param file from {config_url}. "
            f"Cannot proceed with Pygustus without it.\nError: {e}"
        )
        sys.exit(1)
