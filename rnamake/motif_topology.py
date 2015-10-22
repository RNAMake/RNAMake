
import motif_graph
import motif_tree

def graph_to_tree(mg, start=None):
        #if start is None:
            #leafs = mg.leafs()
            #if len(leafs) == 0:
            #    start = mg.graph.oldest_node()

        seen = []
        open = [ start ]

        mt = motif_tree.MotifTree()
        ss = mg.secondary_structure()
        while len(open) > 0:
            current = open.pop(0)
            ss_m = ss.motif(current.data.id)
            print ss_m

            #mt.add_motif(start.data)