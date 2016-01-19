import glob
import subprocess

unittests = glob.glob("*unittest")
for test in unittests:
    subprocess.call("./"+test, shell=True)