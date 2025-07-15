#!/usr/bin/env python3
"""Simple parallel jobs manager using ProcessPoolExecutor."""

import subprocess
import os
import concurrent.futures
import multiprocessing
from toga_modules.common import to_log


class ParallelJobsManager:
    """
    Simple class for managing parallel jobs using ProcessPoolExecutor.
    """

    def __init__(self):
        """Initialize the manager."""
        self.return_code = None
        self.output_data = []

    def execute_jobs(self, joblist_path, manager_data, label, max_workers=None, wait=True):
        """
        Execute jobs in parallel using ProcessPoolExecutor.

        :param joblist_path: Path to the joblist file.
        :param manager_data: Data from the manager class.
        :param label: Label for the run.
        :param max_workers: Maximum number of workers
        :param wait: Boolean -> controls whether run blocking or not
        """
        # Read the joblist to get commands
        commands = []
        with open(joblist_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    commands.append(line)

        to_log(f"ProcessPoolExecutor: executing {len(commands)} jobs")
        
        # Get number of workers or use default
        if max_workers is None:
            max_workers = min(multiprocessing.cpu_count(), len(commands))
        to_log(f"ProcessPoolExecutor: using {max_workers} workers")

        # Create log directory
        log_dir = manager_data["logs_dir"]
        os.makedirs(log_dir, exist_ok=True)
        
        if wait:
            self._execute_jobs(commands, manager_data, label, max_workers)
        else:
            # For non-blocking execution, we'd need to implement this
            # For now, just execute synchronously
            self._execute_jobs(commands, manager_data, label, max_workers)

    def _execute_jobs(self, commands, manager_data, label, max_workers):
        """Execute the jobs using ProcessPoolExecutor."""
        try:
            # Execute all commands in parallel using subprocess
            with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
                # Submit all jobs
                futures = [
                    executor.submit(self._run_single_command, cmd, i, manager_data, label) 
                    for i, cmd in enumerate(commands)
                ]
                
                # Collect results
                results = []
                for i, future in enumerate(concurrent.futures.as_completed(futures)):
                    try:
                        result = future.result()
                        results.append(result)
                        if i % 10 == 0:  # Log progress every 10 completed jobs
                            to_log(f"ProcessPoolExecutor: completed {i+1}/{len(commands)} jobs")
                    except Exception as exc:
                        to_log(f"ProcessPoolExecutor: job failed with exception: {exc}")
                        results.append(None)
                
                self.return_code = 0
                to_log(f"ProcessPoolExecutor: all {len(commands)} jobs completed")
                
        except Exception as e:
            to_log(f"ProcessPoolExecutor: execution failed: {e}")
            self.return_code = 1

    def _run_single_command(self, cmd, job_id, manager_data, label):
        """Run a single command and capture its output."""
        try:
            # Parse the command to extract actual python command and output redirection
            parts = cmd.split(' > ')
            python_cmd = parts[0].strip()
            output_file = parts[1].strip() if len(parts) > 1 else None
            
            # Execute the command
            result = subprocess.run(
                python_cmd,
                shell=True,
                capture_output=True,
                text=True,
                cwd=manager_data.get("project_path", os.getcwd())
            )
            
            # Write output to the specified file if provided
            if output_file and result.stdout:
                with open(output_file, 'w') as f:
                    f.write(result.stdout)
            
            # Log any errors
            if result.stderr:
                log_file = os.path.join(manager_data["logs_dir"], f"{label}_{job_id}.err")
                with open(log_file, 'w') as f:
                    f.write(result.stderr)
            
            return result.returncode == 0
            
        except Exception as e:
            to_log(f"ProcessPoolExecutor: error running command {cmd}: {e}")
            return False

    def check_status(self):
        """Check if the jobs are done."""
        return self.return_code

    def terminate_process(self):
        """Terminate associated process - no-op for ProcessPoolExecutor."""
        pass
