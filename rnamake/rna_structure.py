import motif_type
import x3dna
import structure
import util
import chain_closure
import exceptions
import user_warnings
import residue
import basic_io
import secondary_structure
import basepair
import settings

import primitives.rna_structure
import primitives.basepair
from primitives.rna_structure import ends_from_basepairs, assign_end_id
from primitives.rna_structure import end_id_to_seq_and_db

import os
import numpy as np

class RNAStructure(primitives.rna_structure.RNAStructure):
    """
    Complete container for representing a RNA. Contains both the 3D structure
    information but also includes basepair objects to represent the pairs
    between residues. RNAStructure is rarely called directly but serves as an
    abstract class for both Motif and Pose so for example use please see those
    classes.

    :param struct: The 3D coordinate information
    :param basepairs: basepairing information for each residue pair
    :param ends: the basepairs at the end of chains. These define connection
        points to other RNAStructures
    :param name: name of RNAStructure
    :param path: location of where RNAStructure where 3D information was loaded
        from either the pdb or directory.
    :param mtype: The type of a motif, only needed for Motif.motif objects
    :param score: The score generated by motif_scorer.MotifScorer

    :type struct: structure.Structure
    :type basepairs: list of basepair.Basepairs
    :type ends: list of basepair.Basepairs
    :type name: str
    :type path: str
    :type mtype: motif_type
    :type score: float

    :attributes:

    `structure`: structure.Structure
        structure containing residue and chain information for
        The 3D coordinate information
    `basepairs`: list of basepair.Basepairs
        Basepairs between residues
    `ends`: list of baspir.Basepairs
        Basepair ends where RNA structures can be connected
    `mtype`: motif_type
        The type of a motif, only needed for Motif.motif objects
    `name`: str
        the name of the RNAStructure
    `path`: str
        location of where RNAStructure originated from, this is just
        a place holder for converting from rna_structure.RNAStructure
    `score` : float
        the score generated by motif_scorer.MotifScorer, estimates secondary
        structure stability
    `end_ids`: list of strs
        strings indenifying the secondary structure and sequence
        in the perspective of a given basepair end
    `beads` : list of residue.Bead objects
        keeps the beads required for steric clash calulations

    :examples:

    ..  code-block:: python

        # load test structure
        >>> from rnamake.unittests import instances
        >>> r_struct = instances.rna_structure()

    """

    __slots__ = [
        "_structure",
        "_basepairs",
        "_ends",
        "_name",
        "_score",
        "_protein_beads",
        "_end_ids",
        "_block_end_add",
        "_dot_bracket"
    ]

    def __init__(self, structure, basepairs, ends, end_ids, name, dot_bracket,
                 block_end_add=-1, protein_beads=None):
        self._structure       = structure
        self._basepairs       = basepairs
        self._ends            = ends
        self._end_ids         = end_ids
        self._name            = name
        self._block_end_add   = block_end_add
        self._dot_bracket     = dot_bracket
        self._protein_beads   = protein_beads

        if self._protein_beads is None:
            self._protein_beads = []

        for c in self._structure.get_chains():
            for r in c:
                found = 0
                for i, end in enumerate(ends):
                    if (end.res1_uuid == r.uuid or end.res2_uuid == r.uuid) and \
                       i != self.block_end_add:
                        found = 1
                        break
                if not found:
                    r.build_beads()

    def __iter__(self):
        return self._structure.__iter__()

    @classmethod
    def from_str(cls, s, rts):
        spl = s.split("&")
        name = spl[0]
        block_end_add = int(spl[1])
        struc = structure.Structure.from_str(spl[2], rts)
        bp_strs = spl[3].split("@")
        bps = []
        for bp_str in bp_strs[:-1]:
            bps.append(bp_from_str(struc, bp_str))
        end_strs = spl[4].split("@")
        ends = []
        for end_str in end_strs[:-1]:
            ends.append(bp_from_str(struc, end_str))
        end_ids = spl[5].split()
        bead_strs = spl[6].split(";")
        protein_beads = []
        for bead_str in bead_strs[:-1]:
            protein_beads.append(residue.Bead.from_str(bead_str))

        return cls(struc, bps, ends, end_ids, name, block_end_add, spl[7],
                   protein_beads)

    @classmethod
    def copy(cls, rs, new_uuid=0):
        s = structure.Structure.copy(rs._structure, new_uuid)
        basepairs = []
        ends = []
        protein_beads = [ residue.Bead.copy(b) for b in rs._protein_beads ]

        for bp in rs._basepairs:
            if new_uuid:
                bp_res = rs.get_bp_res(bp)
                res1 = s.get_residue(num=bp_res[0].num, chain_id=bp_res[0].chain_id,
                                     i_code=bp_res[0].i_code)
                res2 = s.get_residue(num=bp_res[1].num, chain_id=bp_res[1].chain_id,
                                     i_code=bp_res[1].i_code)
                bp_new = basepair.Basepair.copy_with_new_uuids(bp, res1.uuid, res2.uuid)
                basepairs.append(bp_new)
            else:
                basepairs.append(basepair.Basepair.copy(bp))

        for end in rs._ends:
            if new_uuid:
                bp_res = rs.get_bp_res(end)
                res1 = s.get_residue(num=bp_res[0].num, chain_id=bp_res[0].chain_id,
                                     i_code=bp_res[0].i_code)
                res2 = s.get_residue(num=bp_res[1].num, chain_id=bp_res[1].chain_id,
                                     i_code=bp_res[1].i_code)
                bp = basepair.Basepair.copy_with_new_uuids(end, res1.uuid, res2.uuid)
                ends.append(bp)
            else:
                ends.append(basepair.Basepair.copy(end))

        return cls(s, basepairs, ends, rs._end_ids, rs._name, rs._block_end_add,
                   rs._dot_bracket, protein_beads)

    def to_str(self):
        """
        stringifies rna structure object
        """
        s  = self._name + "&" + str(self._block_end_add) + "&"
        s += self._structure.to_str() + "&"
        for bp in self._basepairs:
            res1, res2 = self.get_bp_res(bp)
            s += bp.to_str() + ";"
            s += str(res1.num) + "|" + res1.chain_id + "|" + res1.i_code + ";"
            s += str(res2.num) + "|" + res2.chain_id + "|" + res2.i_code + "@"
        s += "&"
        for end in self._ends:
            res1, res2 = self.get_bp_res(end)
            s += end.to_str() + ";"
            s += str(res1.num) + "|" + res1.chain_id + "|" + res1.i_code + ";"
            s += str(res2.num) + "|" + res2.chain_id + "|" + res2.i_code + "@"
        s += "&"
        for end_id in self._end_ids:
            s += end_id + " "
        s += "&"
        s += basic_io.beads_to_str(self._protein_beads)
        s += "&"
        s += self._dot_bracket
        s += "&"
        return s

    def to_pdb_str(self, renumber=-1, close_chain=0):
        """
        creates a PDB string formatted verision of this Structure object.

        :param renumber: what should the first residue be numbered. -1 is
            to NOT renumber, Default=-1.
        :param close_chain: fixes the phosphate backbone, takes a while, so
            default is 0 not to run

        :type renumber: int
        :type close_chain: int

        :return: str
        """
        if close_chain:
            for c in self.iter_chains():
                chain_closure.close_chain(c)

        return self._structure.to_pdb_str(renumber)

    def to_pdb(self, fname="motif.pdb", renumber=-1, close_chain=0):
        """
        write structure to pdb file

        :param fname: name of the file of the pdb file you want to write to
        :param renumber: what should the first residue be numbered. -1 is
            to NOT renumber, Default=-1.
        :param close_chain: fixes the phosphate backbone, takes a while, so
            default is 0 not to run

        :type fname: str
        :type renumber: int
        :type close_chain: int

        :return: None

        """
        #if close_chain:
        #    for c in self._chains():
        #        chain_closure.close_chain(c)

        return self._structure.to_pdb(fname, renumber)

    def sequence(self):
        """
        wrapper for :func:`rnamake.secondary_structure.Structure.sequence`
        """

        seq = ""
        for i, c in enumerate(self._structure.get_chains()):
            if i != 0:
                seq += "&"
            for r in c:
                seq += r.name
        return seq

    def move(self, p):
        self._structure.move(p)
        for bp in self._basepairs:
            bp.move(p)
        for end in self._ends:
            end.move(p)

    def transform(self, t):
        self._structure.transform(t)
        for bp in self._basepairs:
            bp.transform(t)
        for end in self._ends:
            end.transform(t)

    def get_secondary_structure(self):
        db_chains = self._dot_bracket.split("&")
        res = []
        for i, c in enumerate(self._structure.get_chains()):
            for j, r in enumerate(c):
                s_r = secondary_structure.Residue(r.name, db_chains[i][j], r.num,
                                                  r.chain_id, r.i_code, r.uuid)
                res.append(s_r)
        s = secondary_structure.Structure(res, self._structure._chain_cuts)
        bps = []
        for bp in self._basepairs:
            s_bp = secondary_structure.Basepair(bp.res1_uuid, bp.res2_uuid, bp.name, bp.uuid)
            bps.append(s_bp)
        ends = []
        for end in self._ends:
            s_end = secondary_structure.Basepair(end.res1_uuid, end.res2_uuid, end.name, end.uuid)
            ends.append(s_end)

        return secondary_structure.RNAStructure(s, bps, ends, self._end_ids[::])

    def steric_clash(self, rs, clash_radius=settings.CLASH_RADIUS):
        for r1 in self:
            for b1 in r1.iter_beads():
                for r2 in rs:
                    for b2 in r2.iter_beads():
                        if b1.btype == residue.BeadType.PHOS or \
                                        b2.btype == residue.BeadType.PHOS:
                            continue
                        dist = util.distance(b1.center, b2.center)
                        if dist < clash_radius:
                            return 1
        return 0

    @property
    def block_end_add(self):
        return self._block_end_add

    @property
    def dot_bracket(self):
        """
        secondary structure for rna_structure
         """

        return self._dot_bracket



# TODO should probably move this somewhere else?
class ChainEndPairMap(object):
    """
    A simple class container for expressing the relationship between chain
    objects and end objects. An end will always have two residues one at the
    end of a 5' end of a chain and the other at the 3' end of a chain. This
    class organizes the two chains which may be different chains to match the
    two residues in an end. This is primarily used when mergering chains
    together to get the final sequence of a construct.

    :param chain1: the chain that res1 is at the 5' end at in the basepair that this
        object is mapping for
    :param chain2: the chain that res1 is at the 3' end at in the basepair that this
        object is mapping for

    :type chain1: chain.Chain
    :type chain2: chain.Chain

    :attributes:

    `chain1`: chain.Chain
        the chain that res1 is at the 5' end at in the basepair that this
        object is mapping for
    `chain2`: chain.Chain
     the chain that res1 is at the 3' end at in the basepair that this
        object is mapping for

    """

    def __init__(self, chain1=None, chain2=None):
        self.p5_chain, self.p3_chain = chain1, chain2

    def is_hairpin(self):
        """
        A quick check to see if both chains are the same. This would mean
        the chain is a hairpin since both chain ends form a basepair

        :returns: whether both chains are the same
        :rtype: int

        """

        if self.p5_chain == self.p3_chain:
            return 1
        else:
            return 0

    def chains(self):
        """
        returns both chains into an list for easy iteration

        :returns: both chains in a list
        :rtype: list of chain.Chains
        """

        return [self.p5_chain, self.p3_chain]


def rna_structure_from_pdb(pdb_path, rts):
    s = structure.structure_from_pdb(pdb_path, rts)
    bps = basepairs_from_x3dna(pdb_path, s)
    ends = ends_from_basepairs(s, bps)
    end_ids = []

    for end in ends:
        bps.remove(end)

    for end in ends:
        end_id = assign_end_id(s, bps, ends, end)
        end_ids.append(end_id)

    name = util.filename(pdb_path)[:-4]
    score = 0
    seq, dot_bracket = end_id_to_seq_and_db(end_ids[0])
    rna_struc = RNAStructure(s, bps, ends, end_ids, name, dot_bracket=dot_bracket,
                             block_end_add=0)
    return rna_struc


def _calc_center(res):
    center = np.array([0.0,0.0,0.0])
    count = 0
    for r in res:
        for a in r:
            if a is None:
                continue
            center += a.coords
            count += 1
    center /= count
    return center


def bp_from_str(struc, s):
        bp_spl = s.split(";")
        r2_info = bp_spl.pop().split("|")
        r1_info = bp_spl.pop().split("|")
        bp_str = ";".join(bp_spl)
        res1 = struc.get_residue(int(r1_info[0]), r1_info[1], r1_info[2])
        res2 = struc.get_residue(int(r2_info[0]), r2_info[1], r2_info[2])
        return basepair.Basepair.from_str(bp_str, res1.uuid, res2.uuid)


def basepairs_from_x3dna(path, s):
    """
    gets x3dna data on basepairing information and then interwines it
    with the structural information stored in structure for simpler
    retreival of data

    :param path: path to the pdb file
    :param structure: the structure with the same residues that will appear
        in the x3dna output

    :type path: str
    :type structure: structure.Structure

    :return: gets all the basepairs that x3dna finds and returns them as
        basepair.Basepair
    :rtype: list of basepair.Basepair objects

    """
    x3dna_parser = x3dna.X3dna()
    x_basepairs = x3dna_parser.get_basepairs(path)
    basepairs = []
    for xbp in x_basepairs:
        res1 = s.get_residue(num=xbp.res1.num,
                             chain_id=xbp.res1.chain_id,
                             i_code=xbp.res1.i_code)

        res2 = s.get_residue(num=xbp.res2.num,
                             chain_id=xbp.res2.chain_id,
                             i_code=xbp.res2.i_code)

        if res1 is None or res2 is None:
            not_found = 0
            if res1 is None:
                not_found = xbp.res1
            else:
                not_found = xbp.res2
            user_warnings.RNAStructureWarning(
                "cannot find residues in basepair: " + str(not_found) + " "
                "this residue should NOT be a normal nucleotide etc: A,G,C,U\n")
            continue

        try:
            res1.get_atom("C1'")
        except exceptions.ResidueException:
            user_warnings.RNAStructureWarning(
                str(res1) + " has no C1' residue cannot have it in a basepair "
                "is required for alignment\n")
            continue

        try:
            res2.get_atom("C1'")
        except exceptions.ResidueException:
            user_warnings.RNAStructureWarning(
                str(res2) + " has no C1' residue cannot have it in a basepair "
                "is required for alignment\n")
            continue


        center = _calc_center([res1, res2])
        name = primitives.basepair.calc_bp_name([res1, res2])
        bp_type = primitives.basepair.BasepairType.NC

        bp_str = res1.name+res2.name
        wc = "GC,CG,AU,UA".split(",")

        if bp_str in wc and xbp.bp_type == x3dna.X3dnaBPType.cWUW:
            bp_type = primitives.basepair.BasepairType.WC
        elif bp_str == "GU" or bp_str == "UG" and xbp.bp_type == x3dna.X3dnaBPType.cWUW:
            bp_type = primitives.basepair.BasepairType.GU

        bp = basepair.Basepair(res1.uuid, res2.uuid, xbp.r, center,
                               [res1.get_coords("C1'"), res2.get_coords("C1'")],
                               name, bp_type=bp_type, x3dna_bp_type=xbp.bp_type)
        basepairs.append(bp)

    """if os.path.isfile("ref_frames.dat"):
        os.remove("ref_frames.dat")

    if os.path.isfile(name + "_dssr.out"):
        os.remove(name + "_dssr.out")"""

    return basepairs


def get_chain_end_map(chains, end):
    """
    builds :class:`ChainEndPairMap` instances from a target end and a pool
    of chains from a structure.

    :param chains: all chains that could contain residues in the target end
    :param end:

    :type chains: list of chain.Chain objects
    :type end: basepair.Basepair

    :return: a ChainEndPairMap defining the 5' and 3' chains of a given end
    :rtype: ChainEndPairMap
    """

    chain_map = ChainEndPairMap()

    for c in chains:
        for r in end.residues():
            if c.first().uuid == r.uuid and chain_map.p5_chain is None:
                chain_map.p5_chain = c
            elif c.first().uuid == r.uuid:
                raise exceptions.RNAStructureException(
                    "cannot build chain map two residues are assigned to 5' "
                    "chain")

            if c.last().uuid == r.uuid and chain_map.p3_chain is None:
                chain_map.p3_chain = c
            elif c.last().uuid == r.uuid:
                raise exceptions.RNAStructureException(
                    "cannot build chain map two residues are assigned to 3' "
                    "chain")

    if chain_map.p5_chain is None or chain_map.p3_chain is None:
        raise exceptions.RNAStructureException(
            "did not build map properly, both chains are not found")

    return chain_map



























