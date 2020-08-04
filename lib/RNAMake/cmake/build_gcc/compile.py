import argparse
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-c',  help='build preset',
                                    required=True)

    args = parser.parse_args()
    return args


args = parse_args()
subprocess.call("cmake -G Ninja -DCMAKE_CC_COMPILER=gcc-mp-4.9 -DCMAKE_CXX_COMPILER="+args.c, shell=True)
subprocess.call("ninja", shell=True)
