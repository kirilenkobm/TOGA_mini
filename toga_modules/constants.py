"""Class holding all project-wide constants."""
import os

__author__ = "Bogdan M. Kirilenko"
__github__ = "https://github.com/kirilenkobm"

class Constants:
    LOCATION = os.path.dirname(__file__)

    U12_FILE_COLS = 3
    U12_AD_FIELD = {"A", "D"}
    ISOFORMS_FILE_COLS = 2
    NF_DIR_NAME = "nextflow_logs"
    NEXTFLOW = "nextflow"

    PARA_STRATEGIES = ["nextflow", "para", "custom"]  # TODO: add snakemake

    TEMP_CHAIN_CLASS = "temp_chain_trans_class"
    MODULES_DIR = "toga_modules"
    C_LIB_DIR = "util_c/lib"
    C_BIN_DIR = "util_c/bin"
    RUNNING = "RUNNING"
    CRASHED = "CRASHED"
    TEMP = "temp"
