import basic_io
import pdb_parser
import chain

class Structure(object):
    """
	Stores 3D structure information from a pdb file. Stores all chains,
    residues and atoms objects. Implementation is designed to be extremely
    lightweight and capable of performing fast transformations.

    """

    def __init__(self,pdb=None):
        self.chains = []

        if pdb:
            residues = pdb_parser.parse(pdb)
            self._build_chains(residues)

    def _build_chains(self, residues):
        """
		takes all residues and puts into the correct order in chains checking
        for physical connection between O5' and P atoms between residues

		:param residues: residue objects that belong in this structure
		:type residues: List of Residue objects
		"""

        self.chains = []
        #sort residues so check residues for connection quicker as the next on
        #in the array will be closest to it by number
        residues.sort(key=lambda x: x.num)

        while 1:
            current = None
            #find next 5' end, all chains go from 5' to 3'
            for i,r in enumerate(residues):
                five_prime_end = 1
                for j,r2 in enumerate(residues):
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
            #extend chain until 3' end
            while current != None:
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
                    self.chains.append(chain.Chain(current_chain_res))
                    current = None





