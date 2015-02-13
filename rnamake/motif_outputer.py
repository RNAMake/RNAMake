import motif_tree

class MotifOutputer(object):
    def __init__(self):
        self.pdb_str = ""
        self.i = 1
        self.mt = motif_tree.MotifTree(sterics=0)

    def add_motif(self, m, reorient=1):
        if reorient:
            try:
                self.mt.add_motif(m)
            except:
                return

            if len(self.mt.nodes) == 1:
                return

            self.pdb_str += "MODEL "+str(self.i) + "\n"
            self.pdb_str += mt.nodes[1].motif.to_pdb_str() + "\n"
            self.pdb_str += "ENDMDL\n"
        else:
            self.pdb_str += "MODEL "+str(self.i) + "\n"
            self.pdb_str += m.to_pdb_str() + "\n"
            self.pdb_str += "ENDMDL\n"

        self.i += 1

    def to_pdb(self, fname="motifs.pdb"):
        f = open(fname, "w")
        f.write(self.pdb_str)
        f.close()
