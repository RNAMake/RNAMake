import itertools
import x3dna
import structure
import basepair
import util

class Motif(object):
    """
	The basic unit of this project stores the 3D coordinates of a RNA Motif
	as well as the 3DNA parameters such as reference frame and origin for
	each basepair
    """

    def __init__(self, mdir=None, pdb=None):
        self.structure = structure.Structure()
        if mdir is not None and pdb is not None:
            raise ValueError("cannot initiate a Motif with both a mdir and" +\
                             "a pdb")

        if mdir is None and pdb is None:
            self.basepairs = []
            self.mdir = ""
            self.name = ""
            return

        # supplied a motif directory that already contains ref_frames.dat and
        # dssr output file
        if mdir:
            filename = util.filename(mdir)
            self.mdir = mdir
            self.name = filename
            self.structure = structure.Structure(mdir + "/" + filename + ".pdb")

        # supplied only a pdb, have to create the ref_frames and dssr output
        # locally
        if pdb:
            mdir = util.base_dir(pdb)
            filename = util.filename(pdb)
            self.mdir = mdir
            self.name = filename[:-4]
            self.structure = structure.Structure(pdb)

        self.basepairs = self._setup_basepairs()
        self.setup_basepair_ends()

    def _setup_basepairs(self):
        """
        gets x3dna data on basepairing information and then interwines it
        with the structural information stored in structure for simpler
        retreival of data
        """
        x3dna_parser = x3dna.X3dna()
        x_basepairs = x3dna_parser.get_basepairs(self.mdir + "/" + self.name)
        basepairs = []
        for xbp in x_basepairs:
            res1 = self.structure.get_residue(num=xbp.res1.num,
                                              chain_id=xbp.res1.chain_id,
                                              i_code=xbp.res1.i_code)

            res2 = self.structure.get_residue(num=xbp.res2.num,
                                              chain_id=xbp.res2.chain_id,
                                              i_code=xbp.res2.i_code)

            if res1 is None or res2 is None:
                raise ValueError("cannot find residues in basepair")

            bp = basepair.Basepair(res1, res2, xbp.r, xbp.bp_type)
            basepairs.append(bp)

        return basepairs

    def setup_basepair_ends(self):
        # TODO revisit this code, is it really necessary?
        # get chain ends

        chain_ends = []
        for c in self.structure.chains:
            chain_ends.append(c.first())
            if len(c) > 1:
                chain_ends.append(c.last())

        seen_res_bps = { res : [] for res in chain_ends }
        for bp in self.basepairs:
            for ce1 in chain_ends:
                for ce2 in chain_ends:
                    if bp.res1 == ce1 and bp.res2 == ce2:
                        seen_res_bps[ce1].append(bp)
                        seen_res_bps[ce2].append(bp)

        lists = filter(lambda l: len(l) > 0, seen_res_bps.values())
        combos = itertools.product(*lists)
        best, best_count = [], 0
        for c in combos:
            bps, res = [], []
            for bp in c:
                if bp in bps:
                    continue
                fail = 0
                for r in bp.residues():
                    if r in res:
                        fail = 1
                        break
                    res.append(r)
                if fail:
                    continue
                bps.append(bps)
                if len(res) > best_count:
                    best_count = len(res)
                    best = bps

        # for reproducability
        best.sort(key = lambda x : x.name())
        self.ends = best
        return best

    def get_basepair(self, res1=None, res2=None, uuid1=None, uuid2=None):
        found = []
        for bp in self.basepairs:
            if res1 is not None and (res1 != bp.res1 and res1 != bp.res2):
                continue
            if res2 is not None and (res2 != bp.res1 and res2 != bp.res2):
                continue
            if uuid1 is not None and \
               (uuid1 != bp.res1.uuid and uuid1 != bp.res2.uuid):
                continue
            if uuid2 is not None and \
               (uuid2 != bp.res1.uuid and uuid2 != bp.res2.uuid):
                continue
            found.append(bp)
        return found


