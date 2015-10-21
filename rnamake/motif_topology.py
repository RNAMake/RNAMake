
import motif_graph
import motif_tree

def graph_to_tree(mg, start=None):
        if start is None:
            leafs = mg.leafs()
            if len(leafs) == 0:
                start = mg.graph.oldest_node()

        seen = []
        open = [ start ]

        mt = motif_tree.MotifTree()
        while len(open) > 0:
            pass