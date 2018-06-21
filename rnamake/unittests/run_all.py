import glob
import subprocess

unittests = glob.glob("*unittest.py")

for i, unit in enumerate(unittests):
    print unit
    subprocess.call("python " + unit, shell=True)
    #if i > 20:
    #    break