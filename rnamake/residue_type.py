import glob
from . import settings

alt_names = {
    "O1P": "OP1",
    "O2P": "OP2"
}


class ResidueType(object):

    def __init__(self, name, atom_map):
        self.atom_map = atom_map
        self.name = name
        self.alt_names = [name[0], "r"+name[0], "D"+name[0]]


class ResidueTypeSet(object):

    __slots__ = ["residue_types"]

    def __init__(self):
        self.residue_types = []

        path = settings.RESOURCES_PATH + "/residue_types"
        type_files = glob.glob(path + "/*.rtype")
        for tf in type_files:
            name = self._get_rtype_name(tf)
            atom_map = self._get_atom_map_from_file(tf)
            rtype = ResidueType(name, atom_map)
            self.residue_types.append(rtype)

    def _get_rtype_name(self, type_file):
        name_spl = type_file.split("/")
        type_file_name = name_spl[-1]
        return type_file_name[:-6]

    def _get_atom_map_from_file(self, type_file):
        f = open(type_file)
        line = f.readline()
        f.close()
        atom_names = line.split()
        atom_map = {}
        for i, name in enumerate(atom_names):
            atom_map[name] = i
        return atom_map

    def get_rtype_by_resname(self, resname):
        for restype in self.residue_types:
            if resname == restype.name:
                return restype
            if resname in restype.alt_names:
                return restype
        return None


def get_rtype(name):
    return rtypes.get_rtype_by_resname(name)

rtypes = ResidueTypeSet()

# add alt names will at some point be phased out
gua = get_rtype("GUA")
ade = get_rtype("ADE")
ura = get_rtype("URA")
cyt = get_rtype("CYT")

gua.alt_names.extend("MIA GDP GTP M2G 1MG 7MG G7M QUO I YG".split())
ade.alt_names.extend("A23 3DA 1MA 12A AET 2MA".split())
ura.alt_names.extend("PSU H2U 5MU 4SU 5BU 5MC U3H 2MU 70U BRU DT".split())
cyt.alt_names.extend("CBR CCC".split())
