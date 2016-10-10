#!/usr/bin/env python

import argparse
import pandas as pd

from rnamake import resource_manager as rm
from rnamake import util, motif, motif_graph

def parse_args():
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-f', help='path to pdb', required=True)
    args = parser.parse_args()
    return args

args = parse_args()

f = open(args.f)
lines = f.readlines()
f.close()

mg = motif_graph.MotifGraph(mg_str=lines[0])
mg.to_pdb("test.pdb", renumber=1)