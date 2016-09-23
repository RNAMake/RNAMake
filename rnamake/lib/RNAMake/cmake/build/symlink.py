#!/usr/bin/python

import sys
import os

if os.path.exists(sys.argv[2]):
    os.remove(sys.argv[2])

os.symlink(os.path.abspath(sys.argv[1]), sys.argv[2])

