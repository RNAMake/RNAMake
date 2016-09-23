#!/usr/bin/env python

import glob
import subprocess

unittests = glob.glob("../../bin/unittests/*")
for test in unittests:
    subprocess.call(test+" > out", shell=True)
    f = open("out")
    lines = f.readlines()
    f.close()

    # everything passed
    if len(lines) == 3:
        print test, lines[1],
        continue
    else:
        print test, ": FAIL"
        subprocess.call(test, shell=True)


