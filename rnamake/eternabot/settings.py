import os

file_path = os.path.realpath(__file__)
spl = file_path.split("/")
base_dir = "/".join(spl[:-1])

