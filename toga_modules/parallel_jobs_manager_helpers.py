"""Helper function related to parallelisation."""
import time
import os
from toga_modules.common import to_log

__author__ = "Bogdan M. Kirilenko"

ITER_DURATION = 60  # CESAR jobs check interval
NF_DIR_NAME = "nextflow_logs"


def get_nextflow_dir(proj_location, nf_dir_arg):
    """Define nextflow directory."""
    if nf_dir_arg is None:
        default_dir = os.path.join(proj_location, NF_DIR_NAME)
        os.mkdir(default_dir) if not os.path.isdir(default_dir) else None
        return default_dir
    else:
        os.mkdir(nf_dir_arg) if not os.path.isdir(nf_dir_arg) else None
        return nf_dir_arg
