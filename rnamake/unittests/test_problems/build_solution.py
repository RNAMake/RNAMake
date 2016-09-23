import argparse

from rnamake import motif_graph, motif_tree

def parse_args():
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-mg', help='path to motif graph', required=True)
    parser.add_argument('-f', help='path to solution file', required=True)
    args = parser.parse_args()
    return args

args = parse_args()

f = open(args.mg)
lines = f.readlines()
f.close()

f = open(args.f)
sol_lines = f.readlines()
f.close()

mg = motif_graph.MotifGraph(mg_str=lines[0])
start_spl = lines[1].split()
end_spl = lines[2].split()

spl = sol_lines[1].split("\t")
top_str = spl[-1]
mt = motif_tree.motif_tree_from_topology_str(top_str)

mg.add_motif_tree(mt, parent_index=int(start_spl[1]), parent_end_name=start_spl[0])
mg.add_connection(int(end_spl[1]), mg.last_node().index)

f = open("sol.mg", "w")
f.write(mg.to_str() + "\n")
f.close()