import glob
import settings

alt_names = {
    "O1P": "OP1",
    "O2P": "OP2"
}

class ResidueType(object):

    """
    Simple class to hold the topology of the residue, currrently only holds
    which atoms are included in the atomtype but might be later expanded to
    include bonds and charges. This class should not be initiated by itself,
    initiation occurs in ResidueTypeSet

    :param name: residue name
    :param atom_map: the position of where each atom should in a residue by
        name

    :type name: str
    :type atom_map: dict

    :attributes:
    `name` : str
        Residue name
    `atom_map` : dict
        The position of where each atom should in a residue by name
    `alt_name` : list
        Other names the residue can go by, ex. G is also GUA
    """

    __slots__ = ["name", "atom_map", "alt_names"]

    def __init__(self, name, atom_map):
        self.atom_map = atom_map
        self.name = name
        self.alt_names = [name[0], "r"+name[0], "D"+name[0]]

    def __repr__(self):
        return "<ResidueType(name='%s')>" % (self.name)

    def get_correct_atom_name(self, a):
        if a.name in alt_names:
            return [a.name, alt_names[a.name]]
        else:
            return None


class ResidueTypeSet(object):

    """
    Holds all the ResidueType objects, for initiation of new residues. Do not
    initiate a instantance of ResidueTypeSet, if you want a new ResidueType do

    .. code-block:: python

        >>>import rnamake.residue_type
        >>>rnamake.residue_type.get_rtype("GUA")
        <ResidueType(name='GUA')>

    :attributes:
    `residue_types` : list of ResidueTypes
        Contains all residue types that are acceptable in rnamake
    """

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
        """
        extract name from file type_file name.
        :param type_file: file path of residue type file
        :type type_file: str
        """
        name_spl = type_file.split("/")
        type_file_name = name_spl[-1]
        return type_file_name[:-6]

    def _get_atom_map_from_file(self, type_file):
        """
        extracts the atom position for each atom for the Residue object atom
        array for easy indexing
        :param type_file: file path of residue type file
        :type type_file: str
        """
        f = open(type_file)
        line = f.readline()
        f.close()
        atom_names = line.split()
        atom_map = {}
        for i, name in enumerate(atom_names):
            atom_map[name] = i
        return atom_map

    def get_rtype_by_resname(self, resname):
        """
        get the ResidueType object for a given residue by name, this method
        should not be called directly! Should call
        rnamake.residue_type.get_rtype for simplicity
        :param resname: name of residue that you want the ResidueType for
        :type resname: str

        .. code-block:: python
            #get guanine residue type
            >>>rtype =  ResidueTypeSet()
            >>>rtype.get_rtype_by_resname("GUA")
            <ResidueType(name='GUA')>
        """
        for restype in self.residue_types:
            if resname == restype.name:
                return restype
            if resname in restype.alt_names:
                return restype
        return None


def get_rtype(resname):
    """
    Get a reference to a ResidueType by name. This is the only way you should
    get a ResidueType reference

    :param resname: the name of the residue you want the ResidueType
    :type resname: str

    .. code-block:: python

        >>>rnamake.residue_type.get_rtype("GUA")
        <ResidueType(name='GUA')>

    """
    return rtypes.get_rtype_by_resname(resname)

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
