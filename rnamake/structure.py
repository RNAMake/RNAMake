import pdb_parser
import chain
import numpy as np

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

    def __init__(self, pdb=None):
        self.chains = []

        if pdb:
            residues = pdb_parser.parse(pdb)
            self._build_chains(residues)

        self._cache_coords()

    def renumber(self):
        for i, r in enumerate(self.residues()):
            r.num = i+1
            r.chain_id = "A"

    def __repr__(self):
        return "<Structure( Chains: %s>" % (len(self.chains))

    def _build_chains(self, residues):
        """
        takes all residues and puts into the correct order in chains checking
        for physical connection between O5' and P atoms between residues

        :param residues: residue objects that belong in this structure
        :type residues: List of Residue objects
        """

        self.chains = []
        # sort residues so check residues for connection quicker as the next on
        # in the array will be closest to it by number
        residues.sort(key=lambda x: x.num)

        while True:
            current = None
            # find next 5' end, all chains go from 5' to 3'
            for i, r in enumerate(residues):
                five_prime_end = 1
                for j, r2 in enumerate(residues):
                    if r.connected_to(r2) == -1:
                        five_prime_end = 0
                        break
                if five_prime_end:
                    current = r
                    break
            if not current:
                break
            residues.remove(current)
            current_chain_res = []
            # extend chain until 3' end
            while current is not None:
                current_chain_res.append(current)
                found = 0
                for r in residues:
                    if current.connected_to(r) == 1:
                        current = r
                        found = 1
                        break
                if found:
                    residues.remove(current)
                else:
                    # no more residues to add, make chain object
                    self.chains.append(chain.Chain(current_chain_res))
                    current = None

    def _cache_coords(self):
        """
        Stores the atomic coordinates into an array, so can restore atomic
        cordinates after they have been transformed
        """
        atoms = self.atoms()
        self.coords = [a.coords for a in atoms]

    def get_beads(self, excluded_res=[]):
        """
        generates 3-bead model residue beads for all residues in current structure.

        :param excluded_res: List of residue objects that are not to be included in clash checks for collisions
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
            print num,chain_id,i_code
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
        cstruct = Structure()
        for c in self.chains:
            cc = c.copy()
            cstruct.chains.append(cc)

        cstruct._cache_coords()

        return cstruct

    def restore_coords(self):
        self._update_coords(self.coords)

    def _update_coords(self, new_coords):
        for i,a in enumerate(self.atoms()):
            a.coords = new_coords[i]

    def transform(self, t):
        transformed_coords = np.dot(self.coords, t.rotation().T) + \
                             t.translation()

        self._update_coords(transformed_coords)

    def move(self, p):
        for a in self.atoms():
            a.coords += p



