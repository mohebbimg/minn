import subprocess

def run_command():
    command = [
        "python3.10",
        "../scripts/run_trained_MINN_on_testset_file.py",
        "../tests/sample_file_with_pairs_of_microRNA_cts.csv"
    ]
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        print("Command output:\n", result.stdout)
    except subprocess.CalledProcessError as e:
        print("Error occurred while running the command:")
        print(e.stderr)

# Run the function
run_command()
