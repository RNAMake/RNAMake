from rnamake import util

import subprocess

class QSUBJob(object):
    def __init__(self, path, job_str, walltime="1:00:00"):
        self.path = path

        dir = util.base_dir(path)
        s = "#!/bin/bash\n#PBS -o /dev/null\n#PBS -e /dev/null\n#PBS -m n\n#PBS -M nobody@stanford.edu\n#PBS -l walltime=%s\n\n" % (walltime)
        s += "source ~/.bashrc\n"
        s += "cd " + dir + "\n"
        s += job_str

        f = open(path, "w")
        f.write(s)
        f.close()

    def submit(self):
        subprocess.call("qsub " + self.path, shell=True)




