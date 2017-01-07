import uuid
import abc
import exceptions

import util

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


class Transformable(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def move(self, p):
        pass

    @abc.abstractmethod
    def transform(self, t):
        pass


class Residue(BaseStructureObject):
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
        return self._uuid != other._uuid

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


class Chain(BaseStructureObject):
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


class Structure(BaseStructureObject):
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

class Basepair(object):
    __slots__= [
        "_uuid"
    ]

    def __init__(self, bp_uuid=None):
        self._uuid = bp_uuid
        if self._uuid is None:
            self._uuid = uuid.uuid1()

    def __eq__(self, other):
        return self._uuid == other._uuid

    def __ne__(self, other):
        return self._uuid != self._uuid



class RNAStructure(BaseStructureObject):
    __metaclass__ = abc.ABCMeta

    __slots__ = [
        "_structure",
        "_basepairs",
        "_ends",
        "_name",
        "_end_ids"
    ]

    def __init__(self, structure, basepairs, ends, end_ids, name=None):
        self._structure      = structure
        self._basepairs      = basepairs
        self._ends           = ends
        self._name           = name
        self._end_ids        = end_ids

        if self._name is None:
            self._name = ""


    def iter_basepairs(self):
        return self._basepairs.__iter__()

    def iter_ends(self):
        return self._ends.__iter__()

    # wrappers from structure
    def iter_res(self):
        return self._structure.iter_res()

    def iter_chain(self):
        return self._structure.__iter__()

    def get_residue(self, num=None, chain_id=None, i_code=None, uuid=None):
        """
        wrapper for :func:`rnamake.structure.Structure.get_residue`
        """

        return self._structure.get_residue(num=num, chain_id=chain_id,
                                         i_code=i_code, uuid=uuid)

    def num_res(self):
        return self._structure.num_residues()

    def num_chains(self):
        return len(self._structure)

    def get_basepairs(self, bp_uuid=None, res1=None, res2=None, uuid1=None,
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
        for bp in self._basepairs:
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
        for bp in self._basepairs:
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
            if len(found) == 1:
                raise exceptions.RNAStructureException(
                    "more than one basepair has been found with these search criteria")
            found.append(bp)

        if len(found) == 0:
            return None

        return found[0]

    def get_end(self, index=None, end_name=None, end_id=None):
        if index is not None:
            return self._ends[index]

    def get_end_id(self, i):
        return self._end_ids[i]

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

    def get_bp_res(self, bp):
        return [self.get_residue(uuid=bp.res1_uuid), self.get_residue(uuid=bp.res2_uuid) ]

    def num_basepairs(self):
        return len(self._basepairs)

    def num_ends(self):
        return len(self._ends)


class Motif(RNAStructure):
    def __init__(self, struct=None, basepairs=None, ends=None):
        self.structure = struct
        self.basepairs = basepairs
        self.ends = ends


def calc_bp_name(res):
    res1, res2 = res

    res1_name = res1.chain_id+str(res1.num)+str(res1.i_code)
    res2_name = res2.chain_id+str(res2.num)+str(res2.i_code)

    if res1.chain_id < res2.chain_id:
        return res1_name+"-"+res2_name
    if res1.chain_id > res2.chain_id:
        return res2_name+"-"+res1_name

    if res1.num < res2.num:
        return res1_name+"-"+res2_name
    else:
        return res2_name+"-"+res1_name


def ends_from_basepairs(s, bps, check_type=1):
    """
    find basepairs that are composed of two residues who are at the 5' or 3'
    end of their chains. These are elements of alignment where two basepairs
    can be aligned together to build a larger struture.

    :param s: holds all the residues and chains for a structure
    :param bps: All the basepairs extracted with residues in structure,
        generated from :func:`basepairs_from_x3dna`

    :type s: structure.Structure
    :type bps: basepair.Basepair

    :return: the basepairs composed of chain end residues.
    :rtype: list of basepair.Basepairs
    """

    chain_ends_uuids = []
    for c in s:
        chain_ends_uuids.append(c.first().uuid)
        if len(c) > 1:
            chain_ends_uuids.append(c.last().uuid)

    ends = []
    for bp in bps:
        if check_type:
            if bp.bp_type != "cW-W":
                continue
        if not (util.gu_bp(bp, s) or util.wc_bp(bp, s)):
            continue

        if bp.res1_uuid in chain_ends_uuids and bp.res2_uuid in chain_ends_uuids:
            ends.append(bp)

    return ends


def get_res_basepair(basepairs, r):
    for bp in basepairs:
        if bp.res1_uuid == r.uuid or bp.res2_uuid == r.uuid:
            return bp
    return None


def assign_end_id(s, basepairs, end):
    """
    generate a new end_id based on the secondary structure instance in the
    perspective of the supplied end. An end id is a composition of both the
    sequence and secondary structure in a single string.

    Two GC pairs in a row would be: GG_LL_CC_RR. Sequence followed by secondary
    structure with L being left bracket, R being right bracket and U being dot.

    :param ss: secondary structure instance either RNAStructure,Motif or Pose
    :param end: secondary structure basepair that you want the end id to be in
    :return:
    """

    open_chains = []
    for c in s:
        if c.first().uuid == end.res1_uuid or c.first().uuid == end.res2_uuid:
            open_chains.append(c)
            break

    if len(open_chains) == 0:
        raise exceptions.SecondaryStructureException(
            "could not find chain to start with")

    seen_res = {}
    seen_bp = {}
    seen_chains = { open_chains[0] : 1 }
    saved_bp = None
    structure = ""
    seq = ""
    bounds = [0, 0]
    ss_chains = []
    count = 0
    while len(open_chains) > 0:
        c = open_chains.pop(0)
        for r in c:
            count += 1
            dot_bracket = "."
            bp = get_res_basepair(basepairs, r)
            if bp is not None:
                flag =0
                if bp.bp_type != "cW-W":
                    flag = 1
                if not (util.gu_bp(bp, s) or util.wc_bp(bp, s)):
                    flag = 1
                if flag:
                    bp = None
            saved_bp = None
            if bp is not None:
                saved_bp = bp
                partner_res_uuid = bp.partner(r.uuid)
                partner_res = s.get_residue(uuid=partner_res_uuid)
                if   bp not in seen_bp and r not in seen_res and \
                                partner_res not in seen_res:
                    seen_res[r] = 1
                    dot_bracket = "("
                elif partner_res in seen_res:
                    if seen_res[partner_res] > 1:
                        dot_bracket = "."
                    else:
                        dot_bracket = ")"
                        seen_res[r] = 1
                        seen_res[partner_res] += 1

            structure += dot_bracket
            seq += r.name

            if saved_bp is not None:
                seen_bp[saved_bp] = 1

        bounds[1] = count
        ss_chains.append([seq, structure])
        structure = ""
        seq = ""


        best_score = -1

        for c in s:
            if c in seen_chains:
                continue
            score = 0
            for r in c:
                bp = get_res_basepair(basepairs, r)
                if bp is not None and bp in seen_bp:
                    score += 1
            if score > best_score:
                best_score = score

        best_chains = []
        for c in s:
            if c in seen_chains:
                continue
            score = 0
            for r in c:
                bp = get_res_basepair(basepairs, r)
                if bp is not None and bp in seen_bp:
                    score += 1
            if score == best_score:
                best_chains.append(c)

        best_chain = None
        best_score = 10000
        for c in s:
            if c in seen_chains:
                continue
            pos = 1000
            for i, r in enumerate(c):
                bp = get_res_basepair(basepairs, r)
                if bp is not None and bp in seen_bp:
                    pos = i
                    break
            if pos < best_score:
                best_score = pos
                best_chain = c

        if best_chain is None:
            break
        seen_chains[best_chain] = 1
        open_chains.append(best_chain)

    ss_id = ""
    for i, chain in enumerate(ss_chains):
        ss_id += chain[0] + "_"
        for e in chain[1]:
            if   e == "(":
                ss_id += "L"
            elif e == ")":
                ss_id += "R"
            elif e == ".":
                ss_id += "U"
            else:
                raise exceptions.SecondaryStructureException(
                    "unexpected symbol in dot bracket notation: " + e)
        if i != len(ss_chains)-1:
            ss_id += "_"
    return ss_id


def end_id_to_seq_and_db(ss_id):
    ss = ""
    seq = ""
    spl = ss_id.split("_")

    for i in range(0, len(spl)-1, 2):
        seq += spl[i]
        for e in spl[i+1]:
            if   e == "L":
                ss += "("
            elif e == "R":
                ss += ")"
            elif e == "U":
                ss += "."
            else:
                raise ValueError("unexpected symbol in ss_id")

        if i != len(spl)-2:
            seq += "&"
            ss += "&"
    return seq, ss
