import x3dna
import structure
import basepair
import util

class Motif(object):
    def __init__(self,mdir=None,pdb=None):
        self.structure = structure.Structure()

        if mdir:
            filename = util.filename(mdir)
            self.mdir = mdir
            self.name = filename
            self.structure = structure.Structure(mdir + "/" + filename + ".pdb")
            self.basepairs = self._setup_basepairs()

    def _setup_basepairs(self):
        x3dna_parser = x3dna.X3dna()
        x_basepairs = x3dna_parser.get_basepairs(self.mdir + "/" + self.name)
        for xbp in x_basepairs:
            pass
