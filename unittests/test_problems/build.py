from rnamake import motif_graph
from rnamake import resource_manager as rm


def build_mini_ttr():
    mg = motif_graph.MotifGraph()
    ttr = rm.manager.get_motif(name="GAAA_tetraloop", end_name="A229-A245")
    h1 = rm.manager.get_motif(name="HELIX.IDEAL.3")
    h2 = rm.manager.get_motif(name="HELIX.IDEAL.3")
    h3 = rm.manager.get_motif(name="HELIX.IDEAL.3")

    mg.add_motif(h1)
    mg.add_motif(ttr)
    mg.add_motif(h2, parent_end_index=1)
    mg.add_motif(h3, parent_index=1, parent_end_index=2)

    start = mg.get_end(pos=3)
    end = mg.get_end(pos=2)

    f = open("mini_ttr.mg", "w")
    f.write(mg.to_str() + "\n")
    f.write(start.name() + " " + str(3) + "\n")
    f.write(end.name() + " " + str(2) + "\n")
    f.close()


build_mini_ttr()