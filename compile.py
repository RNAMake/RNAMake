import os
import subprocess
import shutil
from rnamake import settings, util

path = util.base_dir(os.path.realpath(__file__))

os.chdir(settings.LIB_PATH + "/lib/RNAMake/cmake/build")

if not os.path.isfile("build.ninja") and os.path.isfile("CMakeCache.txt"):
    os.remove("CMakeCache.txt")
    shutil.rmtree("CMakeFiles")

subprocess.call("python make_project.py", shell=True)
subprocess.call("cmake -G Ninja", shell=True)
subprocess.call("ninja", shell=True)