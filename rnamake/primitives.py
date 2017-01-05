import uuid
import abc
import exceptions

class BaseStructureObject(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def copy(cls, c):
        pass

    @abc.abstractmethod
    def from_str(cls, s, rts):
        pass

    @abc.abstractmethod
    def to_str(self):
        pass


class StructureObject(BaseStructureObject):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def move(self, p):
        pass

    @abc.abstractmethod
    def transform(self, t):
        pass


class Residue(StructureObject):
    __metaclass__ = abc.ABCMeta

    __slots__ = [
        "_name",
        "_num",
        "_chain_id",
        "_i_code",
        "_uuid"]

    def __init__(self, name, num, chain_id, i_code=None, r_uuid=None):
        self._name = name
        self._num = num
        self._chain_id = chain_id
        self._i_code = i_code
        self._uuid = r_uuid

        if self._i_code is None:
            self._i_code = ""

        if self._uuid is None:
            self._uuid = uuid.uuid1()

    def __eq__(self, other):
        return self._uuid == other._uuid

    def __ne__(self, other):
        return self._uuid != self._uuid

    @property
    def name(self):
        return self._name

    @property
    def num(self):
        return self._num

    @property
    def chain_id(self):
        return self._chain_id

    @property
    def i_code(self):
        return self._i_code

    @property
    def uuid(self):
        return self._uuid


class Chain(StructureObject):
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


class Structure(StructureObject):
    __metaclass__ = abc.ABCMeta

    __slots__ = ["_chains"]

    def __init__(self, chains=None):
        self._chains = chains
        if self._chains is None:
            self._chains = []

    def __len__(self):
        return len(self._chains)

    def __iter__(self):
        return self._chains.__iter__()

    def iter_res(self):
        for c in self._chains:
            for r in c:
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

        return found[0]

    def chain(self, i):
        return self._chains[i]

    def num_residues(self):
        total = 0
        for c in self._chains:
            total += len(c)
        return total



class Basepair(object):
    def __init__(self, r, d, sugars, bp_uuid=None):
        self.uuid = bp_uuid



class RNAStructure(object):
    def __init__(self):
        pass


class Motif(RNAStructure):
    def __init__(self, struct=None, basepairs=None, ends=None):
        self.structure = struct
        self.basepairs = basepairs
        self.ends = ends