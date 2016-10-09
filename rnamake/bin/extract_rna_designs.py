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

df = pd.read_csv(args.f)
#print df['mg_str'][0]

for i,r in df.iterrows():
    mg = motif_graph.MotifGraph(mg_str=r['mg_str'])
    print len(mg)
