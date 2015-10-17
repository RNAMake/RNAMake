
class Chain(object):

    """
    Store chain information from pdb file. Stores all residues in chain.
    Implementation is designed to be extremely lightweight.

    :param residues: the residues that are to be included in this chain
    :type residues: list of residue objects

    Attributes:
    ----------
    `residues` : List of Residue objects
        The list of residues that belong to this chain will always be in
        5' to 3' order

    """
    __slots__ = ["residues"]

    def __init__(self, residues=[]):
        self.residues = residues

    def __repr__(self):
        if len(self.residues) == 0:
            return "<Chain( First: None\n\t  Last:  None\n\t  Size: 0)>"

        return "<Chain( First: %s\n\t  Last:  %s\n\t  Size: %s)>" %\
            (self.first(), self.last(), len(self.residues))

    def __len__(self):
        return len(self.residues)

    def first(self):
        """
        returns 5' end of chain
        """

        return self.residues[0]

    def last(self):
        """
        returns 3' end of chain
        """
        return self.residues[-1]

    def subchain(self, start=None, end=None, start_res=None, end_res=None):
        """
        Creates a new chain from a subsection of the residues in the current
        chain.

        :param start: start position in residues object list
        :type start: int

        :param end: end position in residues object list
        :type end: int
        """
        if start_res:
            try:
                start = self.residues.index(start_res)
                end = self.residues.index(end_res)
            except:
                return None

            if start > end:
                start, end = end, start
            end += 1

        if start < 0:
            raise ValueError("start cannot be less then 0")

        if end is None:
            end = len(self.residues)

        return Chain(self.residues[start:end])

    def copy(self):
        """
        Returns a deepcopy of the this chain element
        """
        residues = [r.copy() for r in self.residues]
        return Chain(residues)

    def to_str(self):
        s = ""
        for r in self.residues:
            s += r.to_str() + ";"
        return s

    def to_pdb_str(self, acount=1, return_acount=0):
        s = ""
        for r in self.residues:
            r_str, acount = r.to_pdb_str(acount, 1)
            s += r_str

        if return_acount:
            return s, acount
        else:
            return s

    def to_pdb(self, fname="chain.pdb"):
        f = open(fname, "w")
        f.write(self.to_pdb_str())
        f.close()

    def list_res(self):
        s = ""
        for r in self.residues:
            s += r.rtype.name[0] + str(r.num) + " "
        return s


def connect_residues_into_chains(residues):
        """
        takes all residues and puts into the correct order in chains checking
        for physical connection between O5' and P atoms between residues

        :param residues: residue objects that belong in this structure
        :type residues: List of Residue objects
        """

        chains = []
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
                    chains.append(Chain(current_chain_res))
                    current = None

        return chains

