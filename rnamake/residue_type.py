import glob
import os

import settings, util


alt_names = {
    "O1P": "OP1",
    "O2P": "OP2"
}


class SetType(object):
    RNA = 0,
    PROTEIN = 1


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

    __slots__ = [
        "__name",
        "__atom_map",
        "__alt_names",
        "__set_type"]

    def __init__(self, name, atom_map, set_type, extra_alt_names=None):
        self.__atom_map = atom_map
        self.__set_type = set_type
        self.__name = name
        if self.__set_type == SetType.RNA:
            self.__alt_names = [name[0], "r"+name[0], "D"+name[0]]
        else:
            self.__alt_names = []

        if extra_alt_names is not None:
            self.__alt_names.extend(extra_alt_names)

    def __repr__(self):
        return "<ResidueType(name='%s')>" % (self.__name)

    def __len__(self):
        return len(self.__atom_map)

    def is_valid_atom(self, name):
        if name in self.__atom_map:
            return True
        else:
            return False

    def atom_index(self, name):
        if name in self.__atom_map:
            return self.__atom_map[name]
        else:
            return None

    def get_correct_atom_name(self, a):
        if a.get_name() in alt_names:
            return [a.get_name(), alt_names[a.get_name()]]
        else:
            return None

    def is_alt_name(self, name):
        if name in self.__alt_names:
            return True
        else:
            return False

    @property
    def name(self):
        return self.__name

    @property
    def short_name(self):
        return self.__name[0]

    @property
    def set_type(self):
        return self.__set_type



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

    __slots__ = ["__residue_types"]

    def __init__(self):
        extra_alt_names = {
            'GUA' : "MIA GDP GTP M2G 1MG 7MG G7M QUO I YG".split(),
            'ADE' : "A23 3DA 1MA 12A AET 2MA".split(),
            'URA' : "PSU H2U 5MU 4SU 5BU 5MC U3H 2MU 70U BRU DT".split(),
            'CYT' : "CBR CCC"
        }

        self.__residue_types = []

        path = settings.RESOURCES_PATH + "/residue_types"
        files =  glob.glob(path + "/*")
        for f in files:
            if os.path.isfile(f):
                continue
            rtype_files = glob.glob(f + "/*")
            set_type_name = util.filename(f)
            set_type = None
            if set_type_name == "RNA":
                set_type = SetType.RNA
            elif set_type_name == "PROTEIN":
                set_type = SetType.PROTEIN

            for tf in rtype_files:
                name = self.__get_rtype_name(tf)
                alt_names = None
                if name in extra_alt_names:
                    alt_names = extra_alt_names[name]
                atom_map = self.__get_atom_map_from_file(tf)
                rtype = ResidueType(name, atom_map, set_type, alt_names)
                self.__residue_types.append(rtype)

    def __get_rtype_name(self, type_file):
        """
        extract name from file type_file name.
        :param type_file: file path of residue type file
        :type type_file: str
        """
        name_spl = type_file.split("/")
        type_file_name = name_spl[-1]
        return type_file_name[:-6]

    def __get_atom_map_from_file(self, type_file):
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

    def get_type(self, resname):
        """
        get the ResidueType object for a given residue by name, this method
        should not be called directly! Should call
        rnamake.residue_type.get_rtype for simplicity
        :param resname: name of residue that you want the ResidueType for
        :type resname: str

        .. code-block:: python
            #get guanine residue type
            >>>rtype =  ResidueTypeSet()
            >>>rtype.get_type("GUA")
            <ResidueType(name='GUA')>
        """
        for restype in self.__residue_types:
            if resname == restype.name:
                return restype
            if restype.is_alt_name(resname):
                return restype
        return None

