import abc

from rnamake import util, transform, exceptions
import base
import basepair



class RNAStructure(base.BaseStructureObject):
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

    def __iter__(self):
        return self._structure.__iter__()

    def iter_basepairs(self):
        return self._basepairs.__iter__()

    def iter_ends(self):
        return self._ends.__iter__()

    # wrappers from structure
    def iter_chains(self):
        return self._structure.iter_chains()

    def get_residue(self, num=None, chain_id=None, i_code=None, uuid=None, index=None):
        """
        wrapper for :func:`rnamake.structure.Structure.get_residue`
        """

        return self._structure.get_residue(num, chain_id, i_code, uuid, index)

    def get_res_index(self, r):
        return self._structure.get_res_index(r)

    def get_chains(self):
        return self._structure.get_chains()

    def get_chain(self, i):
        return self._structure.get_chain(i)

    def num_res(self):
        return self._structure.num_residues()

    def num_chains(self):
        return len(self._structure)

    def get_basepairs(self, bp_uuid=None, uuid1=None, uuid2=None, name=None):
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

        if  uuid1 is None and uuid2 is None and bp_uuid is None and name is None:
            raise exceptions.RNAStructureException(
                "no arguments specified for get_basepair()")

        found = []
        for bp in self._basepairs:
            if bp_uuid is not None and bp_uuid != bp.uuid:
                continue
            if uuid1 is not None and \
               (uuid1 != bp.res1_uuid and uuid1 != bp.res2_uuid):
                continue
            if uuid2 is not None and \
               (uuid2 != bp.res1_uuid and uuid2 != bp.res2_uuid):
                continue
            if name is not None and name != bp.name:
                continue
            found.append(bp)
        return found

    def get_basepair(self, index=None, bp_uuid=None, uuid1=None, uuid2=None, name=None):
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
        if  uuid1 is None and uuid2 is None and bp_uuid is None \
                and name is None and index is None:
            raise exceptions.RNAStructureException(
                "no arguments specified for get_basepair()")

        found = []
        for i, bp in enumerate(self._basepairs):
            if index is not None and index != i:
                continue
            if bp_uuid is not None and bp_uuid != bp.uuid:
                continue
            if uuid1 is not None and \
               (uuid1 != bp.res1_uuid and uuid1 != bp.res2_uuid):
                continue
            if uuid2 is not None and \
               (uuid2 != bp.res1_uuid and uuid2 != bp.res2_uuid):
                continue
            if name is not None and name != bp.name:
                continue
            if len(found) == 1:
                raise exceptions.RNAStructureException(
                    "more than one basepair has been found with these search criteria")
            found.append(bp)

        if len(found) == 0:
            return None

        return found[0]

    def get_end(self, index=None, bp_uuid=None, uuid1=None, uuid2=None,
                      name=None, end_id=None):

        if  uuid1 is None and uuid2 is None and bp_uuid is None \
                and name is None and index is None and end_id is None:
            raise exceptions.RNAStructureException(
                "no arguments specified for get_basepair()")

        found = []
        for i, bp in enumerate(self._ends):
            if index is not None and index != i:
                continue
            if end_id is not None and self._end_ids[i] != end_id:
                continue
            if bp_uuid is not None and bp_uuid != bp.uuid:
                continue
            if uuid1 is not None and \
               (uuid1 != bp.res1_uuid and uuid1 != bp.res2_uuid):
                continue
            if uuid2 is not None and \
               (uuid2 != bp.res1_uuid and uuid2 != bp.res2_uuid):
                continue
            if name is not None and name != bp.name:
                continue
            if len(found) == 1:
                raise exceptions.RNAStructureException(
                    "more than one basepair has been found with these search criteria")
            found.append(bp)

        if len(found) == 0:
            return None

        return found[0]

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
            bp = self.get_end(name=name)
            if bp is None:
                raise exceptions.RNAStructureException(
                    "cannot find basepair with name "+name)

            return self._ends.index(bp)
        else:
            matching = []
            for i, end_id in enumerate(self._end_ids):
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

    @property
    def name(self):
        return self._name


def ends_from_basepairs(s, bps):
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
    for c in s.get_chains():
        chain_ends_uuids.append(c.first().uuid)
        if len(c) > 1:
            chain_ends_uuids.append(c.last().uuid)

    ends = []
    for bp in bps:
        if bp.bp_type == basepair.BasepairType.NC:
            continue

        if bp.res1_uuid in chain_ends_uuids and bp.res2_uuid in chain_ends_uuids:
            ends.append(bp)

    return ends


def get_res_basepair(basepairs, r):
    for bp in basepairs:
        if bp.res1_uuid == r.uuid or bp.res2_uuid == r.uuid:
            return bp
    return None


def assign_end_id(s, basepairs, ends, end):
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
    chains = s.get_chains()
    for c in chains:
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
            if bp is None:
                bp = get_res_basepair(ends, r)

            if bp is not None:
                flag =0
                if bp.bp_type == basepair.BasepairType.NC:
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
        for c in chains:
            if c in seen_chains:
                continue
            score = 0
            for r in c:
                bp = get_res_basepair(basepairs, r)
                if bp is None:
                    bp = get_res_basepair(ends, r)
                if bp is not None and bp in seen_bp:
                    score += 1
            if score > best_score:
                best_score = score

        best_chains = []
        for c in chains:
            if c in seen_chains:
                continue
            score = 0
            for r in c:
                bp = get_res_basepair(basepairs, r)
                if bp is None:
                    bp = get_res_basepair(ends, r)
                if bp is not None and bp in seen_bp:
                    score += 1
            if score == best_score:
                best_chains.append(c)

        best_chain = None
        best_score = 10000
        for c in chains:
            if c in seen_chains:
                continue
            pos = 1000
            for i, r in enumerate(c):
                bp = get_res_basepair(basepairs, r)
                if bp is None:
                    bp = get_res_basepair(ends, r)
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


def align_rna_structure(ref_bp, motif_end, motif):
    """
    This is the workhorse of the entire suite. Aligns one end of a motif to
    the reference frame and origin of a Basepair object.

    :param ref_bp: the base pair that the motif end is going to align too
    :param motif_end: the motif end basepair to overly with the ref_bp
    :param motif: the motif object that you want to align

    :type ref_bp: Basepair object
    :type motif_end: Basepair object
    :type motif: Motif object
    """

    r1 , r2 = ref_bp.r , motif_end.r
    r = util.unitarize(r1.T.dot(r2))
    trans = -motif_end.d
    t = transform.Transform(r, trans)
    motif.transform(t)
    bp_pos_diff = ref_bp.d - motif_end.d
    motif.move(bp_pos_diff)

    #alignment is by center of basepair, it can be slightly improved by
    #aligning the c1' sugars
    res1_coord, res2_coord = motif_end.sugars
    ref_res1_coord, ref_res2_coord = ref_bp.sugars

    dist1 = util.distance(res1_coord, ref_res1_coord)
    dist2 = util.distance(res2_coord, ref_res1_coord)

    if dist1 < dist2:
        sugar_diff_1 = ref_res1_coord - res1_coord
        sugar_diff_2 = ref_res2_coord - res2_coord
    else:
        sugar_diff_1 = ref_res1_coord - res2_coord
        sugar_diff_2 = ref_res2_coord - res1_coord

    if dist1 < 5 or dist2 < 5:
        motif.move( (sugar_diff_1 + sugar_diff_2) / 2 )