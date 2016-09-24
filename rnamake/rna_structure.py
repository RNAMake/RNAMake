import motif_type
import basepair
import x3dna
import structure
import util
import chain_closure
import exceptions
import user_warnings

import os
import numpy as np


class RNAStructure(object):
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

    def __init__(self, struct=None, basepairs=None, ends=None, name="assembled",
                 path="assembled", mtype=motif_type.UNKNOWN, score=0):
        self.structure = struct
        if self.structure is None:
            self.structure = structure.Structure()
        self.secondary_structure = None
        self.basepairs = basepairs
        if self.basepairs is None:
            self.basepairs = []
        self.mtype = mtype
        self.name = name
        self.path = path
        self.score = score
        self.ends = ends
        self.protein_beads = []
        if self.ends is None:
            self.ends = []
        self.beads = []
        self.end_ids = []

    def copy(self):
        rna_struct = RNAStructure()
        rna_struct.name      = self.name
        rna_struct.path      = self.path
        rna_struct.score     = self.score
        rna_struct.mtype     = self.mtype
        rna_struct.structure = self.structure.copy()
        rna_struct.beads     = [b.copy() for b in self.beads]
        rna_struct.end_ids   = list(self.end_ids)
        rna_struct.protein_beads = [b.copy() for b in self.protein_beads]

        for bp in self.basepairs:
            new_res1 = rna_struct.get_residue(uuid=bp.res1.uuid)
            new_res2 = rna_struct.get_residue(uuid=bp.res2.uuid)
            # hopefully this doesnt happen anymore
            if new_res1 is None or new_res2 is None:
                raise ValueError("could not find a residue during copy")
            new_r = np.copy(bp.bp_state.r)
            new_bp = basepair.Basepair(new_res1, new_res2, new_r, bp.bp_type)
            new_bp.uuid = bp.uuid
            rna_struct.basepairs.append(new_bp)

        for end in self.ends:
            index = self.basepairs.index(end)
            rna_struct.ends.append(rna_struct.basepairs[index])

        return rna_struct

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
            for c in self.chains():
                chain_closure.close_chain(c)

        return self.structure.to_pdb_str(renumber)

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
        if close_chain:
            for c in self.chains():
                chain_closure.close_chain(c)

        return self.structure.to_pdb(fname, renumber)

    def get_residue(self, num=None, chain_id=None, i_code=None, uuid=None):
        """
        wrapper for :func:`rnamake.structure.Structure.get_residue`
        """

        return self.structure.get_residue(num=num, chain_id=chain_id,
                                          i_code=i_code, uuid=uuid)

    def residues(self):
        """
        wrapper for :func:`rnamake.structure.Structure.residues`
        """

        return self.structure.residues()

    def chains(self):
        """
        wrapper for rnamake.structure.Structure.chain
        """

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

        :examples:

        ..  code-block:: python

            # load test structure
            >>> from rnamake.unittests import instances
            >>> r_struct = instances.rna_structure()

            # get the first basepair for testing purposes
            >>> print r_struct.basepairs[0]
            <Basepair(A1-A24)>

            # retrieve the basepair by name
            >>> r_struct.get_basepair(name="A1-A24")
            [<Basepair(A1-A24)>]

            # retrieve it by a residue in the basepair, either by object
            # reference or unique indentifer.
            >>> res1 = r_struct.basepairs[0].res1
            >>> r_struct.get_basepair(res1=res1)
            [<Basepair(A1-A24)>]

            >>> r_struct.get_basepair(uuid1=res1.uuid)
            [<Basepair(A1-A24)>]

            # Using its indentifer is safer
            # as copying RNA structure will yeild different references
            >>> r_struct_copy = r_struct.copy()
            >>> r_struct_copy.get_basepair(res1=res1)
            []

            >>> r_struct_copy.get_basepair(uuid1=res1.uuid)
            [<Basepair(A1-A24)>]
        """

        if res1 is None and res2 is None and uuid1 is None and uuid2 is None \
           and bp_uuid is None and name is None:
            raise exceptions.RNAStructureException(
                "no arguments specified for get_basepair()")

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
        """
        generates 3-bead model residue beads for all residues in current
        rna_structure.

        :param excluded_res: List of residue objects whose beads are not to be
            included. This is generally end residues that would instantly clash
            with residues they are being overlayed onto when performing motif
            aligning
        :param excluded_ends: List of ends where the resiudes be added to
            exclude list. This is just a fast way to remove residues at
            connection sites between RNA structures.

        :type excluded_res: List of residue.Residue objects
        :type excluded_ends: List of  basepair.Basepair objects

        :return: steric beads for this structure
        :rtype:  List of residue.Bead objects

        :examples:

        ..  code-block:: python

            # load test structure
            >>> from rnamake.unittests import instances
            >>> r_struct = instances.rna_structure()

            # standard use of this function for building
            # this removes sterics from the end basepairs allowing them
            # to be overlayed on another basepair to align too
            >>> r_struct.get_beads(r_struct.ends)

        """

        excluded = []
        if excluded_ends:
            for end in excluded_ends:
                excluded.extend(end.residues())

        if excluded_res:
            excluded.extend(excluded_res)

        self.beads = self.structure.get_beads(excluded) + self.protein_beads
        return self.beads

    def get_end_index(self, name=None, id=None):
        """
        gets the internal end index for an end either by its name or id, not
        used very often.

        :param name: name of end from :func:`rnamake.basepair.Basepair.name`
        :param id: corresponding end id

        :type name: str
        :type id: str

        :returns: index of the end in the internal ends list
        :rtype: int

        """

        if name is None and id is None:
            raise exceptions.RNAStructureException(
                "must specify name or id in get_end_index")

        if name is not None:
            bps = self.get_basepair(name=name)
            if len(bps) == 0:
                raise exceptions.RNAStructureException(
                    "cannot find basepair with name "+name)

            end = bps[0]
            return self.ends.index(end)
        else:
            matching = []
            for i, end_id in enumerate(self.end_ids):
                if end_id == id:
                    matching.append(i)
            if len(matching) > 1:
                raise exceptions.RNAStructureException(
                    "more then one end with id "+ id + " in get_end_index")
            return matching[0]

    def sequence(self):
        """
        wrapper for :func:`rnamake.secondary_structure.Structure.sequence`
        """

        return self.secondary_structure.sequence()

    def dot_bracket(self):
        """
        wrapper for :func:`rnamake.secondary_structure.Structure.dot_bracket`
        """

        return self.secondary_structure.dot_bracket()


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


# TODO dont really need name, all deleting can happen in x3nda module ...
def basepairs_from_x3dna(path, name, structure):
    """
    gets x3dna data on basepairing information and then interwines it
    with the structural information stored in structure for simpler
    retreival of data

    :param path: path to the pdb file
    :param name: the name of structure
    :param structure: the structure with the same residues that will appear
        in the x3dna output

    :type path: str
    :type name: str
    :type structure: structure.Structure

    :return: gets all the basepairs that x3dna finds and returns them as
        basepair.Basepair
    :rtype: list of basepair.Basepair objects

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
            if res1.get_atom("C1'") is None:
                continue
        except exceptions.ResidueException:
            user_warnings.RNAStructureWarning(
                str(res1) + " has no C1' residue cannot have it in a basepair "
                "is required for alignment\n")
            continue

        try:
            if res2.get_atom("C1'") is None:
                continue
        except exceptions.ResidueException:
            user_warnings.RNAStructureWarning(
                str(res2) + " has no C1' residue cannot have it in a basepair "
                "is required for alignment\n")
            continue

        bp = basepair.Basepair(res1, res2, xbp.r, xbp.bp_type)
        basepairs.append(bp)

    if os.path.isfile("ref_frames.dat"):
        os.remove("ref_frames.dat")

    if os.path.isfile(name + "_dssr.out"):
        os.remove(name + "_dssr.out")

    return basepairs


def ends_from_basepairs(structure, basepairs):
    """
    find basepairs that are composed of two residues who are at the 5' or 3'
    end of their chains. These are elements of alignment where two basepairs
    can be aligned together to build a larger struture.

    :param structure: holds all the residues and chains for a structure
    :param basepairs: All the basepairs extracted with residues in structure,
        generated from :func:`basepairs_from_x3dna`

    :type structure: structure.Structure
    :type basepairs: basepair.Basepair

    :return: the basepairs composed of chain end residues.
    :rtype: list of basepair.Basepairs
    """

    chain_ends_uuids = []
    for c in structure.chains:
        chain_ends_uuids.append(c.first().uuid)
        if len(c) > 1:
            chain_ends_uuids.append(c.last().uuid)

    ends = []
    for bp in basepairs:
        if bp.bp_type != "cW-W":
            continue
        if not (util.gu_bp(bp) or util.wc_bp(bp)):
            continue

        if bp.res1.uuid in chain_ends_uuids and bp.res2.uuid in chain_ends_uuids:
            ends.append(bp)

    return ends


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



























