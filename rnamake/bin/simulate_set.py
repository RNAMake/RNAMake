#!/usr/bin/env python

import argparse
import pandas as pd
from rnamake.wrappers.ddg_tecto_wrapper import DDGTectoWrapper

valid_sim_args = {name : 1 for name in "cseq,css,fseq,fss,steps,n".split(",") }

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-csv', help='dataframe in csv format', required=True)
    parser.add_argument('-steps', type=int, help='steps in each simulation')
    parser.add_argument('-n', type=int, help='number of times to run each simulation')
    args = parser.parse_args()
    return args


def is_dataframe_valid(df):
    if not {'cseq', 'css', 'fseq', 'fss'}.issubset(df.columns):
        raise ValueError("supplied csv must contain: cseq, css, fseq and fss")


def get_valid_sim_args(arg_dict):
    sim_args = {}
    for name, val in arg_dict.iteritems():
        if val is None:
            continue
        if name not in valid_sim_args:
            continue
        sim_args[name] = val
    return sim_args


def main():
    args = parse_args()
    df = pd.read_csv(args.csv)
    is_dataframe_valid(df) # check to make sure has all required columns

    # get sim args
    cl_sim_args = get_valid_sim_args(vars(args))

    all_col_names = list(df.columns)
    all_col_names.extend(cl_sim_args.keys())
    all_col_names.append("avg_hit_count")
    all_col_names.append("stdev_hit_count")

    w = DDGTectoWrapper()
    f = open("results.csv", "w")
    f.write(",".join(all_col_names)+"\n")
    for i, row in df.iterrows():
        sim_args = get_valid_sim_args(row.to_dict())
        sim_args.update(cl_sim_args)
        w.run(**sim_args)
        r = w.get_results()

        for col in df.columns:
            f.write(str(row[col]) + ",")
        for key in cl_sim_args.keys():
            f.write(str(cl_sim_args[key]) + ",")
        f.write(str(r.avg) + "," + str(r.stdev) + "\n")

    f.close()




if __name__ == "__main__":
    main()
