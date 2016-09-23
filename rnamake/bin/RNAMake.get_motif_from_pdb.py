#!/usr/bin/env python

import argparse

from rnamake import resource_manager as rm
from rnamake import util, motif

def parse_args():
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-pdb', help='path to pdb', required=True)
    parser.add_argument('-name', required=False)
    parser.add_argument('-n', required=False)
    args = parser.parse_args()
    return args

args = parse_args()

fname = util.filename(args.pdb)

name = fname[:-4]
if args.name:
    name = args.name

rm.manager.add_motif(path=args.pdb, name=name, include_protein=1, align=0)
m = rm.manager.get_motif(name=name)

f = open(name + ".motif", "w")
f.write(m.to_str())
f.close()