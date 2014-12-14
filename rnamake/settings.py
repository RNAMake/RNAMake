import os

file_path = os.path.realpath(__file__)
spl = file_path.split("/")
base_dir = "/".join(spl[:-1])

LIB_PATH = base_dir
RESOURCES_PATH = LIB_PATH + "/resources/"
UNITTEST_PATH = LIB_PATH + "/unittests/"
MOTIF_DIRS = RESOURCES_PATH + "motifs/"
