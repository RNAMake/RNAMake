import os
import platform

file_path = os.path.realpath(__file__)
spl = file_path.split("/")
base_dir = "/".join(spl[:-1])

LIB_PATH = base_dir
RESOURCES_PATH = LIB_PATH + "/resources/"
UNITTEST_PATH = LIB_PATH + "/unittests/"
MOTIF_DIRS = RESOURCES_PATH + "motifs/"
ETERNABOT_PATH = LIB_PATH + "/eternabot/"

OS = None
if platform.system() == 'Linux':
    OS = 'linux'
elif platform.system() == 'Darwin':
    OS = 'osx'
else:
    raise SystemError(platform.system() + " is not supported currently")

VIENNA_BIN = RESOURCES_PATH + "vienna/%s/" % (OS)
X3DNA_PATH = RESOURCES_PATH + "x3dna/%s/" % (OS)
PRECOMPUTED_PATH = RESOURCES_PATH + "precomputed/"
CLASH_RADIUS = 2.5
