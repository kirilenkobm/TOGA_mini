"""Utility class for Toga manager."""
import sys
import platform
from datetime import datetime as dt


class TogaUtil:
    @staticmethod
    def generate_project_name():
        """Generate project name automatically."""
        today_and_now = dt.now().strftime("%Y.%m.%d_at_%H:%M:%S")
        project_name = f"TOGA_project_on_{today_and_now}"
        return project_name

    @staticmethod
    def get_shared_lib_extension():
        """Return the appropriate shared library extension for the current platform."""
        system = platform.system().lower()
        if system == "darwin":  # macOS
            return ".dylib"
        elif system == "windows":
            return ".dll"
        else:  # Linux and other Unix-like systems
            return ".so"

    @staticmethod
    def terminate_parallel_processes(jobs_managers):
        for job_manager in jobs_managers:
            job_manager.terminate_process()