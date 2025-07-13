"""Utility class for Toga manager."""
import sys
from datetime import datetime as dt

from toga_modules.common import to_log


class TogaUtil:
    @staticmethod
    def generate_project_name():
        """Generate project name automatically."""
        today_and_now = dt.now().strftime("%Y.%m.%d_at_%H:%M:%S")
        project_name = f"TOGA_project_on_{today_and_now}"
        return project_name

    @staticmethod
    def log_python_version():
        to_log(f"# python interpreter path: {sys.executable}")
        to_log(f"# python interpreter version: {sys.version}")

    @staticmethod
    def terminate_parallel_processes(jobs_managers):
        to_log(f"KeyboardInterrupt: terminating {len(jobs_managers)} running parallel processes")
        for job_manager in jobs_managers:
            job_manager.terminate_process()