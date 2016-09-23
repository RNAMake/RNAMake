import os
import subprocess
from rnamake import settings, util

path = util.base_dir(os.path.realpath(__file__))



os.chdir(settings.LIB_PATH + "/lib/RNAMake/cmake/build")
subprocess.call("python make_project.py", shell=True)
subprocess.call("cmake -G Ninja", shell=True)
subprocess.call("ninja", shell=True)