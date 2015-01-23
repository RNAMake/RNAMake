from setuptools import setup
import os
import rnamake

print "add this to your .bashrc python path"
print "export PYTHONPATH=$PYTHONPATH:"+ rnamake.util.base_dir(os.path.realpath(__file__))


exit()
setup(name='RNAMake',
      packages=['rnamake']
      )

