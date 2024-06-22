import os


def find_binaries_in_path():
    # Get the PATH environment variable
    path = os.environ.get('PATH', '')

    # Split the PATH into individual directories
    path_dirs = path.split(':')

    # Initialize a list to store the binaries
    binaries = []

    # Iterate through each directory in the PATH
    for directory in path_dirs:
        # Check if the directory exists
        if os.path.isdir(directory):
            # List all files in the directory
            for file in os.listdir(directory):
                file_path = os.path.join(directory, file)
                # Check if the file is executable
                if os.path.isfile(file_path) and os.access(file_path, os.X_OK):
                    binaries.append(file)

    # Remove duplicates and sort the list
    binaries = sorted(list(set(binaries)))

    return binaries


def get_mopac():
    binaries = find_binaries_in_path()
    if "mopac" in binaries:
        return "mopac"
    elif "run_mopac7" in binaries:
        return "run_mopac7"
    else:
        print("No MOPAC binary found, fallback to unstable search")
        for binary in binaries:
            if "mopac" in binary.lower():
                last_binary = binary
                if "run" in binary:
                    return binary
    return last_binary


if __name__ == "__main__":
    print(get_mopac())
