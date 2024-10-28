import os
import warnings


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
    """
    Checks if the MOPAC binary is within the path. looks for approximate
    alternatives. E.g., finds run_mopac7 on ubuntu22, but older MOPAC
    versions might parse files differently, thus being incompatible with
    this wrapper
    """
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
    if "last_binary" not in locals():
        last_binary = None
        print("MOPAC binary not found, returning None")
    return last_binary

    try:
        for lib in libs:
            exec(f"{lib}", globals())
    except:
        raise ImportError()


def optional_imports(libs, namespace):
    if isinstance(libs, str):
        libs = [libs]
    if not isinstance(libs, list):
        raise AttributeError("please supply either a list or str to import")

    for lib in libs:
        try:
            exec(lib, namespace)
        except Exception as e:
            warnings.warn(f"Failed to import {lib}\n{e}", ImportWarning)


def xyz_identifier(geometry):
    try:
        geometry = geometry.split("\n")
        if geometry[0].isnumeric():
            geometry = geometry[1:]

        geometry = "\n".join(geometry).strip().split("\n")

        result = all([xyz_line_checker(line) for line in geometry])
        return result, "\n".join([""]+geometry+[""])
    except:
        return False, None


def xyz_line_checker(line):
    line = line.split()
    if len(line) != 4:
        return False
    if not isinstance(line[0], str):
        return False
    return all([float_checker(x) for x in line[1:]])


def float_checker(string):
    try:
        float(string)
        return True
    except:
        return False


def checkOverlap(mol):
    """
    checks mol object vor overlapping atoms
    """
    threshhold = 0.01
    conf = mol.GetConformer(0)
    atoms = mol.GetNumAtoms()
    overlapping_pairs = []

    for i in range(atoms):
        for j in range(i+1, atoms):
            pos_i = conf.GetAtomPosition(i)
            pos_j = conf.GetAtomPosition(j)

            import numpy as np
            dist = np.sqrt(sum((pos_i[k] - pos_j[k])**2 for k in range(3)))
            if dist < threshhold:
                overlapping_pairs.append((i, j))

    return overlapping_pairs


def BlockToXyz(s: str):
    s = s.strip()
    nindex = s.index("\n")
    if s[:nindex].isnumeric():
        s = s[nindex:].strip()
    return s


if __name__ == "__main__":
    print(get_mopac())
    optional_imports("import numpy as np", globals())
    optional_imports(["import matplotlib.pyplot as plt"], globals())
    optional_imports("obviously wrong", globals())
