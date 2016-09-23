import os
import subprocess
import rnamake

print "add this to your .bashrc python path"
print "export PYTHONPATH=$PYTHONPATH:"+ rnamake.util.base_dir(os.path.realpath(__file__))


