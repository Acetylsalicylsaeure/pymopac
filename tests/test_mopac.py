from time import time
from pymopac import MOPAC_PATH
import subprocess
import os


def test_mopac_water():
    ut = str(time()).replace(".", "")
    temp_path = f"/tmp/pymopac_{ut}"
    try:
        os.mkdir(temp_path)
    except Exception as e:
        print("creating temp folder failed", e)
    with open(f"{temp_path}/water.mop", "w") as file:
        water_input = """PM7
test

O          0.95190        0.00764       -0.08578
H          1.91979        0.05300       -0.09570
H          0.67305        0.91422       -0.28406"""
        file.write(water_input)

    result = subprocess.run(
        [MOPAC_PATH, f"{temp_path}/water.mop"], capture_output=True)
    print(result)
    print(result.stdout)
    print(os.listdir(temp_path))
    assert "ended normally" in result.stdout.decode("utf-8").lower()
