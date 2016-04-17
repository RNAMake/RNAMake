import atom
import residue
import residue_type
import exceptions
import user_warnings

def parse(pdb_file):
    """
    very minimalistc pdb parser, currently does not support multiple MODELS
    in NMR structures but works well for what I need at the moment will be
    expanded in future versions. Currently returns an array of Residue object
    this is because it will mostly be called in in Structure object. Also does
    not parse specified connectivity in CONECT statements. If one is interested
    in parsing into a Structure object, please use structure_from_pdb in the
    structure module.

    :param pdb_file: The path to the PDB formatted file you wish to parse
    :type pdb_file: str

    :return: List of residue.Residue objects

    :examples:

    .. code-block:: python

        >>>import rnamake.unittests.files
        >>>residues = parse(rnamake.unittests.files.P4P6_PDB_PATH)
        >>>len(residues)
        157

        #Not yet sorted into chains, residue 106 is the first in this chain
        >>>residues[0]
        <Residue('G250 chain A')>
    """

    try:
        f = open(pdb_file)
        lines = f.readlines()
        f.close()
    except IOError:
        raise exceptions.PDBParserException("cannot parse pdb file: " +
                                            pdb_file + " as it does not exist")

    coordinates = []
    atomnames = []
    resnames = []
    resnums = []
    chainids = []
    icodes = []

    for line_num, line in enumerate(lines):
        startswith = line[0:6]
        if startswith == 'ATOM  ' or startswith == 'HETATM':
            atomname = line[12:16].strip()
            resname = line[17:21].strip()
            chid = line[21]
            alt = line[16]
            try:
                coords = [float(line[30:38]),
                          float(line[38:46]),
                          float(line[46:54])]
            except:
                raise exceptions.PDBParserError('invalid or missing coordinate(s) at '
                                                'line {0}.'.format(line))

            if len(atomname) == 0:
                user_warnings.warn("line " + str(line_num) + ": no atomname detected",
                                   user_warnings.PDBFormatWarning)

            if len(resname) == 0:
                user_warnings.warn("line " + str(line_num) + ": no resname detected",
                                   user_warnings.PDBFormatWarning)

            atomnames.append(atomname)
            resnames.append(resname)
            chainids.append(chid)
            resnums.append(line[22:26])
            icodes.append(line[26])
            coordinates.append(coords)

        # TODO handle multiple models at some point
        elif startswith == 'MODEL':
            raise exceptions.PDBParserException(pdb_file + " contains NMR MODELS " +
                                                "most likely this is not being parsed " +
                                                "properly")

        elif startswith[:3] == 'END':
            break

    residue_atoms = {}
    for i in range(len(atomnames)):
        if resnames[i] == "HOH":
            continue
        key = resnames[i] + " " + resnums[i] + " " + chainids[i] + " " + icodes[i]
        if key not in residue_atoms:
            residue_atoms[key] = []
        already_has = 0
        for a in residue_atoms[key]:
            if a.name == atomnames[i]:
                already_has = 1
                break
        if already_has:
            continue

        residue_atoms[key].append(atom.Atom(atomnames[i],coordinates[i]))

    residues = []
    for key,res_atoms in residue_atoms.iteritems():
        if len(res_atoms) < 6:
            continue
        spl = key.split()
        rtype = residue_type.get_rtype(spl[0])
        if rtype is None:
            user_warnings.warn("restype " + spl[0] + ": is unknown",
                                user_warnings.PDBFormatWarning)
            continue
        icode = ""
        if len(spl) > 3:
            icode = spl[3]
        r = residue.Residue(rtype, spl[0], int(spl[1]), spl[2], icode)
        r.setup_atoms(res_atoms)
        residues.append(r)

    return residues

