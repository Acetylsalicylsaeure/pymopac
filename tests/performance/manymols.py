import subprocess
import sys
import time
import os

os.chdir("./manymols/")


def run_command(command):
    """
    Run a command and return its exit code and execution time in seconds.
    """
    start_time = time.time()  # Record the start time
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        return e.returncode, None  # Return error code if the command fails
    end_time = time.time()  # Record the end time
    execution_time = end_time - start_time  # Calculate execution time
    return 0, execution_time  # Return success and execution time


def main(smi, count):
    print(f"Calculating molecule {smi} for {count} repetitions")

    # Run and time the Python script
    py_command = ["python3", "./launch.py", smi, count]
    py_exit_code, py_time = run_command(py_command)

    if py_exit_code != 0:
        print("Error running the Python script.")
        return

    # Run and time the Shell script
    sh_command = ["bash", "./launch.sh", count]
    sh_exit_code, sh_time = run_command(sh_command)

    if sh_exit_code != 0:
        print("Error running the Shell script.")
        return

    # Calculate time difference
    time_diff = sh_time - py_time
    print(f"\nShell execution is around {time_diff:.2f} seconds faster")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 time_execution.py <SMILES> <COUNT>")
        sys.exit(1)

    smiles_code = sys.argv[1]
    repetition_count = sys.argv[2]

    main(smiles_code, repetition_count)
