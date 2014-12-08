import sys
import logging
import atom
import residue
import residue_type

logging.basicConfig()
logger = logging.getLogger(__name__)

class PDBParserError(Exception):
    pass


def parse(pdb_file):
    """
    very minimalistc pdb parser, currently does not support multiple MODELS
    in NMR structures but works well for what I need at the moment will be
    expanded in future versions. Currently returns an array of Residue object
    this is because it will mostly be called in in Structure object

    .. code-block:: python
        >>>parser = PDBParser
        >>>residues = parser.parse("p4p6.pdb")

    """

    try:
        f = open(pdb_file)
        lines = f.readlines()
        f.close()
    except IOError:
        raise IOError("cannot parse pdb file: " + pdb_file + " as it does not exist")

    coordinates = []
    atomnames = []
    resnames = []
    resnums = []
    chainids = []
    icodes = []

    for line in lines:
        startswith = line[0:6]
        if startswith == 'ATOM  ' or startswith == 'HETATM':
            atomname = line[12:16].strip()
            resname = line[17:21].strip()
            chid = line[21]
            alt = line[16]
            try:
                coords = [ line[30:38], line[38:46], line[46:54] ]
            except:
                raise PDBParserError('invalid or missing coordinate(s) at '
                                         'line {0}.'.format(line))

            atomnames.append(atomname)
            resnames.append(resname)
            chainids.append(chid)
            resnums.append(line[22:26])
            icodes.append(line[26])
            coordinates.append(coords)

        # TODO handle multiple models at some point
        elif startswith == 'ENDMDL' or startswith[:3] == 'END':
             break

    residue_atoms = {}
    for i in range(len(atomnames)):
        if resnames[i] == "HOH":
            continue
        key = resnames[i] + " " + resnums[i] + " " + chainids[i] + " " + icodes[i]
        if key not in residue_atoms:
            residue_atoms[key] = []
        residue_atoms[key].append(atom.Atom(atomnames[i],coordinates[i]))

    residues = []
    for key,res_atoms in residue_atoms.iteritems():
        spl = key.split()
        rtype = residue_type.get_rtype(spl[0])
        if rtype is None:
            logger.warning(spl[0] + " has no residue type, is being skipped")
            continue
        icode = ""
        if len(spl) > 3:
            icode = spl[3]
        r = residue.Residue(rtype, spl[0], int(spl[1]), spl[2], icode)
        r.setup_atoms(res_atoms)
        residues.append(r)

    return residues

