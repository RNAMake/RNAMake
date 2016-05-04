import motif_type
import basepair
import x3dna
import structure
import util
import chain_closure

import os


class RNAStructure(object):
    def __init__(self, struct=None, basepairs=None, ends=None, name="assembled",
                 path="assembled", mtype=motif_type.UNKNOWN, score=0):
        self.structure = struct
        if self.structure is None:
            self.structure = structure.Structure()
        self.basepairs =basepairs
        if self.basepairs is None:
            self.basepairs = []
        self.mtype = mtype
        self.name = name
        self.path = path
        self.score = score
        self.ends = ends
        if self.ends is None:
            self.ends = []
        self.beads = []
        self.end_ids = []

    def to_pdb_str(self, renumber=-1, close_chain=0):
        """
        returns pdb formatted string of motif's structure object
        """
        if close_chain:
            for c in self.chains():
                chain_closure.close_chain(c)

        return self.structure.to_pdb_str(renumber)

    def to_pdb(self, fname="motif.pdb", renumber=-1, close_chain=0):
        """
        writes the current motif's structure to a pdb
        """
        if close_chain:
            for c in self.chains():
                chain_closure.close_chain(c)

        return self.structure.to_pdb(fname, renumber)

    def get_residue(self, num=None, chain_id=None, i_code=None, uuid=None):
        """
        wrapper to self.structure.get_residue()
        """
        return self.structure.get_residue(num=num, chain_id=chain_id,
                                          i_code=i_code, uuid=uuid)

    def residues(self):
        """
        wrapper to self.structure.residues()
        """
        return self.structure.residues()

    def chains(self):
        return self.structure.chains

    def get_basepair(self, bp_uuid=None, res1=None, res2=None, uuid1=None,
                     uuid2=None, name=None):
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
        alt_name = None
        if name:
            name_spl = name.split("-")
            alt_name = name_spl[1] + "-" + name_spl[0]

        found = []
        for bp in self.basepairs:
            if bp_uuid is not None and bp_uuid != bp.uuid:
                continue
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
            if name is not None and name != bp.name():
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

    def get_end_index(self, name=None, id=None):
        if name is None and id is None:
            raise ValueError("must specify name or id in get_end_index")
        if name is not None:
            bps = self.get_basepair(name=name)
            if len(bps) == 0:
                raise ValueError("cannot find basepair with name "+name)
            end = bps[0]
            return self.ends.index(end)
        else:
            matching = []
            for i, end_id in enumerate(self.end_ids):
                if end_id == id:
                    matching.append(i)
            if len(matching) > 1:
                raise ValueError("more then one end with id "+ id + " in get_end_index")
            return matching[0]

    def sequence(self):
        return self.secondary_structure.sequence()

    def dot_bracket(self):
        return self.secondary_structure.dot_bracket()


# TODO should probably move this somewhere else?
class ChainEndPairMap(object):
    def __init__(self, chain1=None, chain2=None):
        self.p5_chain, self.p3_chain = chain1, chain2

    def is_hairpin(self):
        if self.p5_chain == self.p3_chain:
            return 1
        else:
            return 0

    def chains(self):
        return [self.p5_chain, self.p3_chain]


def basepairs_from_x3dna(path, name, structure):
    """
    gets x3dna data on basepairing information and then interwines it
    with the structural information stored in structure for simpler
    retreival of data
    """
    x3dna_parser = x3dna.X3dna()
    x_basepairs = x3dna_parser.get_basepairs(path)
    basepairs = []
    for xbp in x_basepairs:
        res1 = structure.get_residue(num=xbp.res1.num,
                                     chain_id=xbp.res1.chain_id,
                                     i_code=xbp.res1.i_code)

        res2 = structure.get_residue(num=xbp.res2.num,
                                     chain_id=xbp.res2.chain_id,
                                     i_code=xbp.res2.i_code)

        if res1 is None or res2 is None:
            raise ValueError("cannot find residues in basepair")

        bp = basepair.Basepair(res1, res2, xbp.r, xbp.bp_type)
        basepairs.append(bp)

    if os.path.isfile("ref_frames.dat"):
        os.remove("ref_frames.dat")

    if os.path.isfile(name + "_dssr.out"):
        os.remove(name + "_dssr.out")

    return basepairs


def ends_from_basepairs(structure, basepairs):
        chain_ends = []
        for c in structure.chains:
            chain_ends.append(c.first())
            if len(c) > 1:
                chain_ends.append(c.last())

        ends = []
        for bp in basepairs:
            if bp.bp_type != "cW-W":
                continue
            if not (util.gu_bp(bp) or util.wc_bp(bp)):
                continue

            if bp.res1 in chain_ends and bp.res2 in chain_ends:
                ends.append(bp)

        return ends


def get_chain_end_map(chains, end):
    chain_map = ChainEndPairMap()

    for c in chains:
        #for r in c.residues:
        #    print r,
        #print
        for r in end.residues():
            if c.first().uuid == r.uuid and chain_map.p5_chain == None:
                chain_map.p5_chain = c
            elif c.first().uuid == r.uuid:
                raise ValueError("cannot build chain map two residues are assigned"
                                 "to 5' chain")

            if c.last().uuid == r.uuid and chain_map.p3_chain == None:
                chain_map.p3_chain = c
            elif c.last().uuid == r.uuid:
                raise ValueError("cannot build chain map two residues are assigned"
                                 "to 3' chain")

    if chain_map.p5_chain == None or chain_map.p3_chain == None:
        raise ValueError("did not build map properly, both chains are not found")

    return chain_map



























