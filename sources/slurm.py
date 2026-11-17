"""Slurm tasks through submitit, with a record of what was submitted and how it ended."""
import datetime
import getpass
import subprocess
import sys
from pathlib import Path


def submitted(name: str, out: Path, account: str | None, function, *arguments, memory: int = 8, minutes: int = 60) -> list:
    """Runs `function` over `arguments`, a Slurm task each, under the name `name`, and returns the tasks' results. The
    logs go to submitit-logs/ under `out`, with a record of what was submitted and how it ended."""
    import submitit  # only on O2: the group "qc"

    logs = out / "submitit-logs"
    folder, record = logs / name, logs / f"{name}.log"
    if folder.exists() or record.exists() or subprocess.run(["squeue", "-h", "-u", getpass.getuser(), "-n", name], capture_output=True, text=True).stdout.strip():
        sys.exit(f"{name} has run already, or is running: choose another --name")
    logs.mkdir(parents=True, exist_ok=True)
    executor = submitit.AutoExecutor(folder=folder)
    executor.update_parameters(name=name, slurm_partition="short", cpus_per_task=1, mem_gb=memory, timeout_min=minutes, slurm_srun_args=["--export=ALL"], slurm_array_parallelism=20,
                               slurm_additional_parameters={"account": account} if account else {})
    jobs = executor.map_array(function, *arguments)

    def note(text):
        with record.open("a") as log:
            log.write(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} {text}\n")
        print(text, flush=True)

    note(f"submitted {jobs[0].job_id.split('_')[0]}: {len(jobs)} tasks; short, 1 CPU, {memory} GB, {minutes} min")
    submitit.helpers.monitor_jobs(jobs)
    failed = [job.job_id for job in jobs if job.state != "COMPLETED"]
    note(f"finished: {len(jobs) - len(failed)} of {len(jobs)} tasks completed" + (f"; failed: {', '.join(failed)}" if failed else ""))
    return [job.result() for job in jobs]  # raises what a task raised
