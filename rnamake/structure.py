import numpy as np

import pdb_parser
import chain
import util

class StructureException(Exception):
    pass

class Structure(object):
    """
    Stores 3D structure information from a pdb file. Stores all chains,
    residues and atoms objects. Implementation is designed to be extremely
    lightweight and capable of performing fast transformations.

    :param pdb: the path of the pdb to load this structure from, optional
    :type pdb: str

    Attributes
    ----------
    `chains` : list of Chain objects that belong to the current structure
    """

    def __init__(self, chains=[]):
        self.chains = chains
        self.name = "N/A"

    def __repr__(self):
        return """<Structure(name: %s, #chains: %s, #residues: %s, #atoms: %s)>""" %\
               (self.name, len(self.chains), len(self.residues()), len(self.atoms()))

    def renumber(self):
        chains = "A B C D E F G H I J K L M".split()
        j = 1
        for i, c in enumerate(self.chains):
            for r in c.residues:
                r.num = j
                r.chain_id = chains[i]
                j += 1

    def get_beads(self, excluded_res=[]):
        """
        generates 3-bead model residue beads for all residues in current structure.

        :param excluded_res: List of residue objects that are not to be included in clash
         checks for collisions
        :type excluded_res: List of Residue objects
        """

        beads = []
        for r in self.residues():
            if r in excluded_res:
                continue
            beads.extend(r.get_beads())
        return beads

    def get_residue(self, num=None, chain_id=None, i_code=None, uuid=None):
        """
        find a residue based on residue num, chain_id, insert_code and uuid
        will return an error if more then one residue matches search to avoid
        confusion

        :param num: residue number
        :param chain id: what chain the residue belongs to
        :param i_code: the insertation code of the residue
        :param uuid: the unique indentifier that each residue is given

        :type num: int
        :type chain_id: str
        :type i_code: str
        :type uuid: uuid
        """

        found = []
        for c in self.chains:
            for r in c.residues:
                if uuid is not None and uuid != r.uuid:
                    continue
                if num is not None and num != r.num:
                    continue
                if i_code is not None and i_code != r.i_code:
                    continue
                if chain_id is not None and chain_id != r.chain_id:
                    continue
                found.append(r)

        if len(found) > 1:
            self.to_pdb()
            #for f in found:
            #    print f.name, f.num
            #print "num,chain_id,icode,uuid=",num,chain_id,i_code,uuid
            raise ValueError(
                "found multiple residues in get_residue(), narrow " +
                "your search")

        if len(found) == 0:
            return None

        return found[0]

    def residues(self):
        residues = []
        for c in self.chains:
            residues.extend(c.residues)
        return residues

    def atoms(self):
        atoms = []
        for r in self.residues():
            for a in r.atoms:
                if a is None:
                    continue
                atoms.append(a)
        return atoms

    def to_str(self):
        s = ""
        for c in self.chains:
            s += c.to_str() + ":"
        return s

    def to_pdb_str(self):
        acount = 1
        s = ""
        for c in self.chains:
            c_str, acount = c.to_pdb_str(acount, 1)
            s += c_str
            s += "TER\n"
        return s

    def to_pdb(self, fname="structure.pdb"):
        """
        write structure to pdb file

        :param fname: name of the file of the pdb file you want to write to
        :type fname: str

        """
        f = open(fname, "w")
        f.write(self.to_pdb_str())
        f.close()

    def copy(self):
        """
        creates a deep copy of this structure
        """
        chains = []
        for c in self.chains:
            cc = c.copy()
            chains.append(cc)

        return Structure(chains)

    def transform(self, t):
        r_T = t.rotation().T
        for a in self.atoms():
            a.coords = np.dot(a.coords, r_T) + t.translation()

        #transformed_coords = np.dot(self.coords, t.rotation().T) + \
        #                     t.translation()

    def move(self, p):
        for a in self.atoms():
            a.coords += p


def structure_from_pdb(pdb_path):
    residues = pdb_parser.parse(pdb_path)
    chains = chain.connect_residues_into_chains(residues)
    s = Structure(chains)
    s.name = util.filename(pdb_path[:-4])
    return s
