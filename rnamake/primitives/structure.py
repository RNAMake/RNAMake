import uuid
import abc
from rnamake import exceptions

import base

class Structure(base.BaseStructureObject):
    __metaclass__ = abc.ABCMeta

    __slots__ = [
        "_chains",
        "_residues"]

    def __init__(self, chains):
        self._chains = chains
        self._residues = []
        for c in chains:
            for r in c:
                self._residues.append(r)


    def __len__(self):
        return len(self._chains)

    def __iter__(self):
        return self._chains.__iter__()

    def iter_res(self):
        for r in self._residues:
            yield r

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

        if num is None and chain_id is None and i_code is None and uuid is None:
            raise exceptions.StructureException(
                "must specify a parameter to find a residue in get_residue")


        found = []
        for c in self._chains:
            for r in c:
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

    def chain(self, i):
        return self._chains[i]

    def num_residues(self):
        total = 0
        for c in self._chains:
            total += len(c)
        return total
