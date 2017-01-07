import uuid
import abc
from rnamake import exceptions

import base

class Chain(base.BaseStructureObject):
    __metaclass__ = abc.ABCMeta

    __slots__ = ["_residues"]

    def __init__(self, residues=None):
        self._residues = []
        if residues is not None:
            self._residues = residues

    def __len__(self):
        return len(self._residues)

    def __iter__(self):
        return self._residues.__iter__()

    def first(self):
        """
        returns residue at 5' end of chain

        :returns: first residue in chain
        :rtype: residue.Residue

        :examples:

        ..  code-block:: python

            >>> import rnamake.unittests.instances
            >>> c = rnamake.unittests.instances.chain()
            >>> c.first()
            <Residue('G103 chain A')>

        """
        if len(self._residues) == 0:
            raise exceptions.ChainException("cannot call first there are no "
                                            "residues in chain")
        return self._residues[0]

    def last(self):
        """
        returns 3' end of chain
        """
        if len(self._residues) == 0:
            raise exceptions.ChainException("cannot call first there are no "
                                            "residues in chain")
        return self._residues[-1]

    def subchain(self, start=None, end=None, start_res=None, end_res=None):
        """
        Creates a new chain from a subsection of the residues in the current
        chain.

        :param start: start position in residues object list, default:None
        :param end: end position in residues object list, default:None
        :param start_res: The 5' residue of sub chain, default:None
        :param end_res:  The 3' resiude of sub chain, default:None

        :type start: int
        :type end: int

        :return: Chain object

        :examples:

        ..  code-block:: python

            >>> cs = c.subchain(1, 10)
            >>> len(cs)
            9

            >>> cs.first()
            <Residue('A104 chain A')>

            >>> cs2 = c.subchain(start_res=c.residues[10], end_res=c.residues[15])
            >>> len(cs2)
            6

        """

        if start_res is not None and end_res is not None:
            try:
                start = self._residues.index(start_res)
                end = self._residues.index(end_res)
            except:
                raise exceptions.ChainException("supplied start_res and end_res "
                                                "but they are not members of "
                                                "chain")

            if start > end:
                start, end = end, start
            end += 1

        elif start_res is not None and end_res is None:
            raise exceptions.ChainException("supplied start_res but not end_res")

        elif start_res is not None and start is not None:
            raise exceptions.ChainException("cannot supply start and start_res")

        if start < 0:
            raise exceptions.ChainException("start pos cannot be less then 0")

        if end is None:
            end = len(self._residues)


        return self.__class__(self._residues[start:end])

    def contain_res(self, r):
        for res in self._residues:
            if res == r:
                return 1

        return 0

    def residue(self, index):
        return self._residues[index]
