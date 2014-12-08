import residue


class Chain(object):
    """
    Store chain information from pdb file. Stores all residues in chain.
	Implementation is designed to be extremely lightweight.

	Attributes
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

    def subchain(self, start, end):
        """
        Creates a new chain from a subsection of the residues in the current
        chain.

        :param start: start position in residues object list
        :type start: int

        :param end: end position in residues object list
        :type end: int
        """
        if start < 0:
            raise ValueError("start cannot be less then 0")

        return Chain(self.residues[start:end])

    def copy(self):
        """
		Returns a deepcopy of the this chain element
		"""
        residues = [r.copy() for r in self.residues]
        return Chain(residues)


