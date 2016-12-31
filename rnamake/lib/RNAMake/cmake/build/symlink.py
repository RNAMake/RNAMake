#!/usr/bin/python

import sys
import os
import subprocess
import time

if os.path.exists(sys.argv[2]):
    subprocess.call('rm ' + sys.argv[2], shell=True)
    time.sleep(0.50)

os.symlink(os.path.abspath(sys.argv[1]), sys.argv[2])

