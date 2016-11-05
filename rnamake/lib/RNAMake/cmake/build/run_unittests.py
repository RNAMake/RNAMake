#!/usr/bin/env python

import glob
import subprocess
import os

unittests = glob.glob("../../bin/unittests/*")
for test in unittests:

    lines = subprocess.check_output(test, shell=True).split("\n")

    # everything passed
    if len(lines) < 6:
        print test, lines[1]+"\n",
        continue
    elif lines[0:3] == "All":
        print test, lines[1]+"\n",
        continue
    else:
        print test, ": FAIL"
        subprocess.call(test, shell=True)


