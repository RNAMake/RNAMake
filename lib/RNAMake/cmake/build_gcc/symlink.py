#!/usr/bin/python

import sys
import os

if os.path.exists( "../../bin/unittests/"+sys.argv[1]):
    os.remove( "../../bin/unittests/"+sys.argv[1])

os.symlink(os.path.abspath(sys.argv[1]),
           "../../bin/unittests/"+sys.argv[1])

