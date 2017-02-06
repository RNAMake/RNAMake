import uuid
import abc
from rnamake import exceptions

import base

class Structure(base.BaseStructureObject):
    __metaclass__ = abc.ABCMeta

    __slots__ = [
        "_residues",
        "_chain_cuts"]

    def __init__(self, residues, chain_cuts):
        self._residues = residues
        self._chain_cuts = chain_cuts


    def __len__(self):
        return len(self._chain_cuts)

    def __iter__(self):
        return self._residues.__iter__()

    def iter_chains(self):
        for r in self.get_chains():
            yield r

    def get_residue(self, num=None, chain_id=None, i_code=None, uuid=None, index=None):
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

        if num is None and chain_id is None and i_code is None and uuid is None \
           and index is None:
            raise exceptions.StructureException(
                "must specify a parameter to find a residue in get_residue")


        found = []
        for i, r in enumerate(self._residues):
            if index is not None and index != i:
                continue
            if uuid is not None and uuid != r.uuid:
                continue
            if num is not None and num != r.num:
                continue
            if i_code is not None and i_code != r.i_code:
                continue
            if chain_id is not None and chain_id != r.chain_id:
                continue
            found.append(r)

        if len(found) == 0:
            return None

        if len(found) > 1:
            raise exceptions.StructureException(
                "more than one residue was found with this query")

        return found[0]

    def get_res_index(self, res):
        for i, r in enumerate(self._residues):
            if res == r:
                return i
        raise exceptions.StructureException("cannot find res: " + r)

    def get_chains(self):
        raise ValueError("not implmemented")

    def get_chain(self, i):
        return self.get_chains()[i]

    def num_residues(self):
        return len(self._residues)
