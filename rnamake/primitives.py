

class Residue(object):
    def __init__(self, name, num, chain_id, i_code=""):
        self.name = name,
        self.num = num
        self.chain_id = chain_id
        self.i_code = i_code

    def copy(self):
        pass

    def to_str(self):
        pass


class Basepair(object):
    def __init__(self, res1, res2, bp_type="c..."):
        self.res1 = res1
        self.res2 = res2
        self.bp_type = bp_type

    def residues(self):
        return [self.res1, self.res2]

    def name(self):
        str1 = self.res1.chain_id+str(self.res1.num)+str(self.res1.i_code)
        str2 = self.res2.chain_id+str(self.res2.num)+str(self.res2.i_code)

        if str1 < str2:
            return str1+"-"+str2
        else:
            return str2+"-"+str1


class Chain(object):
    def __init__(self, residues=None):
        self.residues = []
        if residues is not None:
            self.residues = residues

    def copy(self):
        pass

    def first(self):
        return self.residues[0]

    def last(self):
        return self.residues[-1]

    def to_str(self):
        pass


class Structure(object):
    def __init__(self, chains=None):
        self.chains = chains
        if self.chains is None:
            self.chains = []
        self.name = "N/A"

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

        if len(found) == 0:
            return None

        return found[0]

    def residues(self):
        residues = []
        for c in self.chains:
            residues.extend(c.residues)
        return residues

    def to_str(self):
        pass

    def copy(self):
        pass


class RNAStructure(object):
    def __init__(self):
        pass


class Motif(RNAStructure):
    def __init__(self, struct=None, basepairs=None, ends=None):
        self.structure = struct
        self.basepairs = basepairs
        self.ends = ends