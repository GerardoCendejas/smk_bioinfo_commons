# This file contains the function definitions to get safe output and log paths
# They return string values (paths) that are in the correspondent safe
# directories (results/ and results/logs/)

# Define base directories
RESULTS_DIR = "results"
LOGS_DIR = "results/logs"
BENCHMARKS_DIR = "results/benchmarks"

# workflow/rules/common.smk

# General rules for getting output and log paths based on the configuration file

def get_output(rule_config_key, ext):
    """
    It formats a secure path for result files based on the configuration.
    Args:
        rule_config_key: The key in the config file for the specific rule (e.g. "windows_phyml")
        I will use the local rule name...
        ext: the extension of the file to be used... Include the dot in the ext (.txt).
    """
    # Extracts the subdir name from config
    subdir = config[rule_config_key]["dir"]
    filename = f"{config[rule_config_key]['prefix']}_{ext}"
    
    # Builds safe path
    return os.path.join(RESULTS_DIR, subdir, filename)

def get_output_filename(rule_config_key, filename_key, ext):
    """
    It formats a secure path for result files based on the configuration.
    Different from get_output, this function allows specifying a different
    filename key in the config file.
    Args:
        rule_config_key: The key in the config file for the specific rule (e.g. "windows_phyml")
        I will use the local rule name...
        filename_key: The key in the config file for the specific filename (e.g. "filename" or "prefix")
        ext: the extension of the file to be used... Include the dot in the ext (.txt).
    """
    # Extracts the subdir name from config
    subdir = config[rule_config_key]["dir"]
    filename = f"{config[rule_config_key][filename_key]}_{ext}"
    
    # Builds safe path
    return os.path.join(RESULTS_DIR, subdir, filename)

def get_log(rule_config_key, wildcard=None):
    """
    It formats a secure path for log files based on the configuration.
    Args:
        rule_config_key: The key in the config file for the specific rule (e.g. "windows_phyml")
        I will use the local rule name, so that results/logs/ follow the same structure as results/...
        wildcard: optional file prefix (no extension) to support wildcarding.
    """
    # Extracts the subdir name from config
    subdir = config[rule_config_key]["dir"]
    prefix = config[rule_config_key]["prefix"]
    filename = f"{prefix}_{wildcard}.log" if wildcard is not None else f"{prefix}.log"

    # Builds safe path
    return os.path.join(LOGS_DIR, subdir, filename)


def get_benchmark(rule_config_key, wildcard=None):
    """
    It formats a secure path for benchmark files based on the configuration.
    Args:
        rule_config_key: The key in the config file for the specific rule (e.g. "windows_phyml")
        wildcard: optional file prefix (no extension) to support wildcarding.
    """
    # Extracts the subdir name from config
    subdir = config[rule_config_key]["dir"]
    prefix = config[rule_config_key]["prefix"]
    filename = f"{prefix}_{wildcard}.tsv" if wildcard is not None else f"{prefix}.tsv"

    # Builds safe path
    return os.path.join(BENCHMARKS_DIR, subdir, filename)


# For envs calling

from pathlib import Path

def _find_commons_root(start: Path) -> Path:
    """Sube hasta encontrar el dir que contiene envs/ (= commons/workflow)."""
    for d in [start, *start.parents]:
        if (d / "envs").is_dir():
            return d
    raise FileNotFoundError(f"No encontré envs/ subiendo desde {start}")

COMMONS_WORKFLOW = _find_commons_root(Path(str(workflow.current_basedir)))
COMMONS_ENVS = COMMONS_WORKFLOW / "envs"
COMMONS_SCRIPTS = COMMONS_WORKFLOW / "scripts"

def get_env(name):
    """Ruta absoluta a un env de commons."""
    p = COMMONS_ENVS / f"{name}.yaml"
    if not p.is_file():
        raise FileNotFoundError(f"Env no existe: {p}")
    return str(p)

def get_script(rel_path):
    """Ruta absoluta a un script de commons (e.g. 'newick_utils/reroot_tree.R')."""
    p = COMMONS_SCRIPTS / rel_path
    if not p.is_file():
        raise FileNotFoundError(f"Script no existe: {p}")
    return str(p)

### Functions that are be useful in some cases

def lines2string(path):
    """
    It takes a path for a file and returns a python list of each line of the file as a string.

    Eg: if you have samples defined in file, you can get them to use as wildcards with this function.
    """
    with open(path, 'r') as f:
        lines = f.readlines()
    return [line.strip() for line in lines]

