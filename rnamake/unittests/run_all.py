import glob
import subprocess

unittests = glob.glob("*unittest.py")

for unit in unittests:
    print unit
    subprocess.call("python " + unit, shell=True)