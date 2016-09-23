#!/usr/bin/env python

import sys
import os
import argparse
import rnamake.prediction.motif_ensemble as motif_ensemble
import rnamake.motif as motif

def parse_args():
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-p', help='path', required=True)
    parser.add_argument('-max', required=False)
    parser.add_argument('-n', required=False)
    args = parser.parse_args()
    return args


args = parse_args()
me = motif_ensemble.MotifEnsemble(args.p, 0 ,0)

max_num = 10
name = "cluster"
if args.n:
    name = args.n

for i, ms in enumerate(me.motif_states):
    m = motif.str_to_motif(ms.mts.build_string)
    m.to_pdb(name+"."+str(i)+".pdb")
    if i > max_num-1:
        break

