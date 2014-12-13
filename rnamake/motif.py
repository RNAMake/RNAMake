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

        self.beads = []
        self.score = 0

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
            self._assign_bp_primes(bp)
            basepairs.append(bp)

        return basepairs

    def _assign_bp_primes(self, bp):
        """
		This is legacy code not sure if I still need it The purpose is to
        determine which residue in a Basepair object is going in the 5' and 3'
        direction. How this is calculated is position from a chain end. If
        residue 1 is positioned in the half the chain its called the 5' end and
        if its in the second half its 3'. If its exactly at the 1/2 point then
        the second residue is checked in the same manner.

		:param bp: an end basepair
		:type bp: Basepair object
		"""
        res1_pos, res2_pos, res1_total, res2_total = 0, 0, None, None
        for c in self.structure.chains:
            for i,r in enumerate(c.residues):
                if bp.res1 == r:
                    res1_pos = i
                    res1_total = len(c.residues)
                elif bp.res2 == r:
                    res2_pos = i
                    res2_total = len(c.residues)

        if  res1_pos > res1_total/2 or res2_pos < res2_total/2:
			bp.res1,bp.res2 = bp.res2,bp.res1

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
        """
        locates a Basepair object based on residue objects or uuids if nothing
        is supplied you will get back all the basepairs in the motif. The way
        to make sure you will only get one basepair back is to supply BOTH
        res1 and res2 OR uuid1 and uuid2, I have left it open like this
        because it is sometimes useful to get all basepairs that a single
        residue is involved

        :param res1: First residue
        :param res2: Second residue
        :param uuid1: First residue uuid
        :param uuid2: Second residue uuid

        :type res1: Residue object
        :type res2: Residue object
        :type uuid1: uuid object
        :type uuid2: uuid object
        """

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

    def get_beads(self, excluded_ends=None, excluded_res=None):
        excluded = []
        if excluded_ends:
            for end in excluded_ends:
                excluded.extend(end.residues())

        if excluded_res:
            excluded.extend(excluded_res)

        self.beads = self.structure.get_beads(excluded)
        return self.beads

    def to_str(self):
        s = self.mdir + "&" + str(self.score) + "&" + \
            self.structure.to_str() + "&"

        #for bp in self.basepairs:


    def to_pdb(self, fname="motif.pdb"):
        pass

