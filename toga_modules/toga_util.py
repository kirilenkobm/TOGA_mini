"""Utility class for Toga manager."""
from datetime import datetime as dt


class TogaUtil:
    @staticmethod
    def generate_project_name():
        """Generate project name automatically."""
        today_and_now = dt.now().strftime("%Y.%m.%d_at_%H:%M:%S")
        project_name = f"TOGA_project_on_{today_and_now}"
        return project_name

    @staticmethod
    def terminate_parallel_processes(jobs_managers):
        for job_manager in jobs_managers:
            job_manager.terminate_process()