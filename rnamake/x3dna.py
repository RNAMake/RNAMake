import subprocess
import settings
import os
import re
import numpy as np
import util
import shutil

import motif_type, exceptions

X3DNA_BIN_PATH = settings.X3DNA_PATH + "bin/"
os.environ['X3DNA'] =  settings.X3DNA_PATH

#- = U
#+ = P
#. = D
class X3dnaBPType(object):
    cmU  = 0  #cm-
    cMUM = 1  #cM-M
    tWPW = 2  #tW+W
    cDPM = 3  #c.+M
    DWPW = 4  #.W+W
    tWUM = 5  #tW-M
    tmUM = 6  #tm-M
    cWPM = 7  #cW+M
    DWUW = 8  #.W-W
    cMPD = 9  #cM+.
    cDUm = 10 #c.-m
    cMPW = 11 #cM+W
    tMPm = 12 #tM+m
    tMUW = 13 #tM-W
    cmUm = 14 #cm-m
    cMUW = 15 #cM-W
    cWUW = 16 #cW-W
    cDUM = 17 #c.-M
    cmPM = 18 #cm+M
    cmUM = 19 #cm-M
    DDDD = 20 #....
    cmUW = 21 #cm-W
    tMUm = 22 #tM-m
    cDUW = 23 #c.-W
    cMPm = 24 #cM+m
    cMUm = 25 #cM-m
    cDDD = 26 #c...
    tWPm = 27 #tW+m
    cDPm = 28 #c.+m
    tmPm = 29 #tm+m
    tWPD = 30 #tW+.
    tmPW = 31 #tm+W
    tDDD = 32 #t...
    cWUD = 33 #cW-.
    cWUM = 34 #cW-M
    tDUW = 35 #t.-W
    tMPM = 36 #tM+M
    tDUM = 37 #t.-M
    cMUD = 38 #cM-.
    cWUm = 39 #cW-m
    tDPm = 40 #t.+m
    tMUD = 41 #tM-.
    cmPW = 42 #cm+W
    cMPM = 43 #cM+M
    cmPD = 44 #cm+.
    cmUD = 45 #cm-.
    cDUD = 46 #c.-.
    cWPW = 47 #cW+W
    tDUD = 48 #t.-.
    tDPW = 49 #t.+W
    tmUm = 50 #tm-m
    cWPD = 51 #cW+.
    tmPD = 52 #tm+.
    tDPD = 53 #t.+.
    cDPD = 54 #c.+.
    tDUm = 55 #t.-m
    tDPM = 56 #t.+M


def get_x3dna_bp_type(name):
    if   name == "cm-":
        return X3dnaBPType.cmU
    elif name == "cM-M":
        return X3dnaBPType.cMUM
    elif name == "tW+W":
        return X3dnaBPType.tWPW
    elif name == "c.+M":
        return X3dnaBPType.cDPM
    elif name == ".W+W":
        return X3dnaBPType.DWPW
    elif name == "tW-M":
        return X3dnaBPType.tWUM
    elif name == "tm-M":
        return X3dnaBPType.tmUM
    elif name == "cW+M":
        return X3dnaBPType.cWPM
    elif name == ".W-W":
        return X3dnaBPType.DWUW
    elif name == "cM+.":
        return X3dnaBPType.cMPD
    elif name == "c.-m":
        return X3dnaBPType.cDUm
    elif name == "cM+W":
        return X3dnaBPType.cMPW
    elif name == "tM+m":
        return X3dnaBPType.tMPm
    elif name == "tM-W":
        return X3dnaBPType.tMUW
    elif name == "cm-m":
        return X3dnaBPType.cmUm
    elif name == "cM-W":
        return X3dnaBPType.cMUW
    elif name == "cW-W":
        return X3dnaBPType.cWUW
    elif name == "c.-M":
        return X3dnaBPType.cDUM
    elif name == "cm+M":
        return X3dnaBPType.cmPM
    elif name == "cm-M":
        return X3dnaBPType.cmUM
    elif name == "....":
        return X3dnaBPType.DDDD
    elif name == "cm-W":
        return X3dnaBPType.cmUW
    elif name == "tM-m":
        return X3dnaBPType.tMUm
    elif name == "c.-W":
        return X3dnaBPType.cDUW
    elif name == "cM+m":
        return X3dnaBPType.cMPm
    elif name == "cM-m":
        return X3dnaBPType.cMUm
    elif name == "c...":
        return X3dnaBPType.cDDD
    elif name == "tW+m":
        return X3dnaBPType.tWPm
    elif name == "c.+m":
        return X3dnaBPType.cDPm
    elif name == "tm+m":
        return X3dnaBPType.tmPm
    elif name == "tW+.":
        return X3dnaBPType.tWPD
    elif name == "tm+W":
        return X3dnaBPType.tmPW
    elif name == "t...":
        return X3dnaBPType.tDDD
    elif name == "cW-.":
        return X3dnaBPType.cWUD
    elif name == "cW-M":
        return X3dnaBPType.cWUM
    elif name == "t.-W":
        return X3dnaBPType.tDUW
    elif name == "tM+M":
        return X3dnaBPType.tMPM
    elif name == "t.-M":
        return X3dnaBPType.tDUM
    elif name == "cM-.":
        return X3dnaBPType.cMUD
    elif name == "cW-m":
        return X3dnaBPType.cWUm
    elif name == "t.+m":
        return X3dnaBPType.tDPm
    elif name == "tM-.":
        return X3dnaBPType.tMUD
    elif name == "cm+W":
        return X3dnaBPType.cmPW
    elif name == "cM+M":
        return X3dnaBPType.cMPM
    elif name == "cm+.":
        return X3dnaBPType.cmPD
    elif name == "cm-.":
        return X3dnaBPType.cmUD
    elif name == "c.-.":
        return X3dnaBPType.cDUD
    elif name == "cW+W":
        return X3dnaBPType.cWPW
    elif name == "t.-.":
        return X3dnaBPType.tDUD
    elif name == "t.+W":
        return X3dnaBPType.tDPW
    elif name == "tm-m":
        return X3dnaBPType.tmUm
    elif name == "cW+.":
        return X3dnaBPType.cWPD
    elif name == "tm+.":
        return X3dnaBPType.tmPD
    elif name == "t.+.":
        return X3dnaBPType.tDPD
    elif name == "c.+.":
        return X3dnaBPType.cDPD
    elif name == "t.-m":
        return X3dnaBPType.tDUm
    elif name == "t.+M":
        return X3dnaBPType.tDPM
    else:
        raise ValueError("cannot get x3dna type with: " + name)


# TODO check to see what happens with RNA without chain ids
class X3dna(object):
    """
    a simple wrapper for interfacing with the x3dna package for determining the
    basepairing in a pdb structure. Please support the author of x3dna directory
    and cite his work and register on his site: http://x3dna.org/

    Data is returned as a list of X3DNA basepair objects which are simplistic
    containers for storing which residues are contained in a basepair and
    what the reference frame and type of the basepair is. These are different
    from basepair.Basepair objects which store the actual atomic information
    of a basepair

    .. code-block:: python

        >>> from rnamake.unittests import files
        >>> from rnamake import x3dna
        >>> x = x3dna.X3dna()
        >>> basepairs = x.get_basepairs(files.P4P6_PDB_PATH)
        >>> len(basepairs)
        79

        >>> basepairs[10]
        <X3DNA Basepair(A131 A192 cW-W)>

    """

    def __init__(self):
        self.find_pair_path = X3DNA_BIN_PATH + "find_pair"
        self.analyze_path   = X3DNA_BIN_PATH + "analyze"
        self.dssr_path      = X3DNA_BIN_PATH + "x3dna-dssr"

        if not os.path.isfile(self.find_pair_path):
            raise exceptions.X3dnaException(
                "%s executable should exist but doesnt! " % (self.find_pair_path))

        if not os.path.isfile(self.analyze_path):
            raise exceptions.X3dnaException(
                "%s executable should exist but doesnt! " % (self.analyze_path))

        if not os.path.isfile(self.dssr_path):
            raise exceptions.X3dnaException(
                "%s executable should exist but doesnt! " % (self.dssr_path))

        self.ref_frames_extra_files = (
            "auxiliary.par,bestpairs.pdb,bp_helical.par,bp_order.dat," + \
            "bp_step.par,cf_7methods.par,col_chains.scr,col_helices.scr," + \
            "hel_regions.pdb,hstacking.pdb,poc_haxis.r3d,stacking.pdb").split(",")

        self.dssr_extra_files = (
            "dssr-2ndstrs.ct,dssr-2ndstrs.dbn,dssr-helices.pdb,"+ \
            "dssr-pairs.pdb,dssr-stems.pdb,hel_regions.pdb,hstacking.pdb,"+ \
            "poc_haxis.r3d,stacking.pdb,dssr-torsions.dat,dssr-Kturns.pdb,"+ \
            "dssr-multiplets.pdb,dssr-hairpins.pdb,dssr-Aminors.pdb").split(",")

    def _remove_files(self, f_list):
        for f in f_list:
            try:
                os.remove(f)
            except:
                pass

    def generate_ref_frame(self, pdb_path):
        """
        runs x3dna/find_pair and x3dna/analyze to get the basepairing
        information of a given pdb and it is updated to a ref_frames.dat
        file which can be parsed by get_basepair_info if necessary

        :param pdb_path: the pdb you want to run without the .pdb extension
        :type pdb_path: str

        .. code-block:: python

            # after running should be a "ref_frames.dat" in current dir
            >>> from rnamake import x3dna
            >>> x = x3dna.X3dna()
            >>> x.generate_ref_frame("test.pdb")

        """
        pdb_name = util.filename(pdb_path)[:-4]
        if not os.path.isfile(pdb_path):
            raise exceptions.X3dnaException(
                pdb_path + " is not found cannot generate ref_frames.dat file")

        result = subprocess.call(
                    self.find_pair_path + " " + pdb_path + " 2> /dev/null " +
                    "stdout | " + self.analyze_path + " stdin >& /dev/null",
                    executable="/bin/bash", shell=True)

        if result != 0:
            raise exceptions.X3dnaException("find_pair did not run correctly")

        if not os.path.isfile("ref_frames.dat"):
            raise exceptions.X3dnaException(
                "generate_ref_frame was ran on %s but no "  % (pdb_path) +
                "ref_frames.dat file was produced something extremely wrong")

        self._remove_files(self.ref_frames_extra_files)
        os.remove(pdb_name+".out")

    def generate_dssr_file(self, pdb_path):
        """
        runs x3dna's dssr which assigns the secondary structure and each of
        the motifs in a RNA structure

        :param pdb_name: the pdb you want to run without the .pdb extension
        :type pdb_name: str

        """
        if not os.path.isfile(pdb_path):
            raise exceptions.X3dnaException(
                pdb_path + " is not found cannot generate dssr file")

        pdb_name = util.filename(pdb_path)[:-4]
        result = subprocess.call(
            self.dssr_path + " -i="+pdb_path+" -o="+pdb_name+ \
            "_dssr.out --non-pair >& /dev/null", shell=True,
            executable="/bin/bash")

        if result != 0:
            raise exceptions.X3dnaException("dssr did not run correctly")

        self._remove_files(self.dssr_extra_files)

    def _parse_ref_frame_file(self, ref_frames_path):
        """
        Helper function do not run outside! parses out basepair information
        from ref_frames.dat which can be generated using generate_ref_frame

        :param ref_frames_path: the path to the ref_frames.dat file
        :type ref_frames_path: str
        """
        f = open(ref_frames_path)
        lines = f.readlines()
        f.close()

        p = re.compile('(\w+)\:\.*(\d+)\_\:\[\.*\w+\]\w+ \- (\w+)\:\.*(\d+)')
        r = []
        d = []
        basepairs = []
        start_bp = 0
        bps = None
        for l in lines:
            if l[0:3] == "...":
                start_bp = 1
                if p.search(l):
                    m = p.search(l)
                    bps = m.groups()
                else:
                    print l
                    raise exceptions.X3dnaException(
                        "did not parse 3DNA ref frame correctly")
                continue
            if start_bp == 0:
                continue
            if start_bp == 1:
                split = re.split("\s+", l)
                d = [ float(x) for x in split[1:4] ]
                start_bp += 1
            elif start_bp > 1:
                split = re.split("\s+", l)
                r.append([ float(x) for x in split[1:4] ])
                start_bp += 1
            #have everything, generate BasePair object
            if start_bp == 5:
                res1 = Residue(int(bps[1]), bps[0])
                res2 = Residue(int(bps[3]), bps[2])
                basepairs.append(Basepair(res1, res2, np.array(r), np.array(d)))
                r = []
                d = []

        return basepairs

    def _get_ref_frame_path(self, path):
        """
        determines were the ref_frames.dat file generated by x3dna is located
        checks where the pdb is where would be the it would be if its a motif
        directory. Also checks local directory and if does not exist in either
        location it creates a new one using generate_ref_frame. This is a
        helper function do not call directly!

        :param pdb_name: the pdb you want to run without the .pdb extension
        :type pdb_name: str

        """
        filename = util.filename(path)

        if os.path.isdir(path):
            ref_frames_path = path + "/ref_frames.dat"
            if not os.path.isfile(ref_frames_path):
                self.generate_ref_frame(path + "/" + filename + ".pdb")

        else:
            base_dir = util.base_dir(path)
            filename = filename[:-4]
            ref_frames_path = base_dir + "/ref_frames.dat"

        if not os.path.isfile(ref_frames_path):
            ref_frames_path = "ref_frames.dat"
            self.generate_ref_frame(path)

        return ref_frames_path

    def _get_dssr_file_path(self, path):
        """
        determines were the pdb_dssr.out file generated by x3dna is located
        checks where the pdb is where would be the it would be if its a motif
        directory. Also checks local directory and if does not exist in either
        location it creates a new one using generate_dssr_file. This is a
        helper function do not call directly!

        :param pdb_name: the pdb you want to run without the .pdb extension
        :type pdb_name: str
        """
        filename = util.filename(path)
        dssr_file_path = None

        if os.path.isdir(path):
            dssr_file_path = path + "/" + filename + "_dssr.out"
            if not os.path.isfile(dssr_file_path):
                self.generate_dssr_file(path + "/" + filename + ".pdb")
                shutil.move(filename + "_dssr.out", path + "/")
        else:
            base_dir = util.base_dir(path)
            filename = filename[:-4]
            dssr_file_path = base_dir + "/" + filename + "_dssr.out"

        if not os.path.isfile(dssr_file_path):
            dssr_file_path = filename + "_dssr.out"
            self.generate_dssr_file(path)

        return dssr_file_path

    def _parse_dssr_res_str(self, res_str):
        """
        helper function to convert the dssr residue indentifier into a
        X3DNA residue object

        :param res_str: dssr residue indentifier
        :type res_str: str

        """
        # TODO reevalute this code, not sure I need all this complexity
        p = re.compile(r"^\d+$")

        spl = re.split("\.",res_str)

        chain = spl[0]
        rnum = spl[1][1:]
        i_code = " "

        if not p.match(rnum):
            r_spl = spl[1]
            i_code = ""
            i_spl = re.split("\^",spl[1])

            if len(i_spl) > 1:
                r_spl = i_spl[0]
                i_code = i_spl[1]

            num = ""
            for e in r_spl[::-1]:
                try:
                    int_e = int(e)
                    num = e + num
                except:
                    break
        try:
            int(rnum)
        except:
            return None


        return Residue(int(rnum), chain, i_code)

    def _divide_dssr_file_into_sections(self, dssr_file_path):
        """
        splits up the dssr file by the title

        :param dssr_file_path: the location of the dssr file
        :type dssr_file_path: str
        """
        f = open(dssr_file_path)
        lines = f.readlines()
        f.close()

        p = re.compile(r"List of \d+\s*(\S+)")

        sections = {}
        section = []
        section_name = ""
        for l in lines:
            if len(section) == 0 and p.match(l):
                m = p.match(l)
                section_name = m.group(1)
                if section_name == 'helix':
                    section_name = 'helices'
                if section_name[-1] != "s":
                    section_name = section_name + "s"
            if l[0] == '*':
                if len(section) != 0 and section_name not in sections:
                    sections[section_name] = section
                section = []
                continue
            section.append(l)

        return sections

    def get_basepairs(self, pdb_path):
        """
        generates and/or collects the information from the ref_frame.dat
        file and dssr output file to determine the basepairing between
        residues and what type each basepair is. Returns a list of basepair
        X3DNA objects.

        :param pdb_path: path to the pdb file
        :type pdb_path: str
        """
        ref_frames_path = self._get_ref_frame_path(pdb_path)
        dssr_file_path = self._get_dssr_file_path(pdb_path)

        basepairs = self._parse_ref_frame_file(ref_frames_path)
        try:
            section = self._divide_dssr_file_into_sections(dssr_file_path)['bases']
        except:
            return []

        for l in section:
            spl = re.split("\s+",l)
            if len(spl) < 6:
                continue
            if not spl[1].isdigit():
                continue
            res1 = self._parse_dssr_res_str(spl[2])
            res2 = self._parse_dssr_res_str(spl[3])
            if res1 is None or res2 is None:
                continue
            bp_type_str = spl[8]
            if len(bp_type_str) < 1:
                bp_type_str = spl[7]

            bp_type = get_x3dna_bp_type(bp_type_str)

            found = 0
            for bp in basepairs:
                if bp.res1 == res1 and bp.res2 == res2:
                    bp.bp_type = bp_type
                    found = 1
                if bp.res2 == res1 and bp.res1 == res2:
                    bp.bp_type = bp_type
                    found = 1

                if found:
                    break

            if not found:
                basepairs.append(Basepair(res1, res2, np.eye(3),
                                          np.array([-1,-1,-1])))
                basepairs[-1].bp_type = bp_type

        return basepairs

    def get_dssr_file_sections(self, pdb_path):
        dssr_file_path = self._get_dssr_file_path(pdb_path)
        return self._divide_dssr_file_into_sections(dssr_file_path)

    def get_motifs(self, pdb_path):
        dssr_file_path = self._get_dssr_file_path(pdb_path)
        types = { 'hairpins'  : motif_type.HAIRPIN,
                  'bulges'    : motif_type.TWOWAY,
                  'internals' : motif_type.TWOWAY,
                  'junctions' : motif_type.NWAY,
                  'non-loops' : motif_type.SSTRAND}

        file_sections = self._divide_dssr_file_into_sections(dssr_file_path)
        all_motifs = []
        for section, t in types.iteritems():
            if section in file_sections:
                motifs = self._parse_dssr_section(file_sections[section], t)
                all_motifs.extend(motifs)
        if 'stems' in file_sections:
            motifs = self._parse_dssr_helix_section(file_sections['stems'])
            all_motifs.extend(motifs)

        return all_motifs

    def _parse_dssr_section(self, section, mtype):
        motifs = []
        seen_res = []
        for l in section:
            spl = l.split()
            try:
                if spl[0][:3] != 'nts':
                    continue
            except:
                continue

            if len(spl) < 3:
                continue
            res = []
            for res_str in spl[2].split(","):
                res_obj = self._parse_dssr_res_str(res_str)
                if res_obj is None:
                    continue
                res.append(res_obj)

            count = 0
            for r in res:
                if r in seen_res:
                    count += 1
            if count == len(res):
                continue
            seen_res.extend(res)
            motifs.append(Motif(res, mtype))

        return motifs

    def _parse_dssr_helix_section(self, section):
        p = re.compile("\d+")
        res = []
        motifs = []
        for l in section:
            spl = l.split()
            try:
                if not p.match(spl[0]):
                    continue
            except:
                continue
            if int(spl[0]) == 1 and len(res) > 0:
                motifs.append(Motif(res, motif_type.HELIX))
                res = []
            r1 = self._parse_dssr_res_str(spl[1])
            r2 = self._parse_dssr_res_str(spl[2])
            if r1 is not None:
                res.append(r1)
            if r2 is not None:
                res.append(r2)

        if len(res) > 0:
            motifs.append(Motif(res, motif_type.HELIX))

        return motifs


class Motif(object):
    """
    Container for motif residues from dssr files produceed by x3dna. Should not
    be instantiated outside x3dna.X3dna class.

    :attributes:
    `residues` : list of Residue objects
        residues that belong to same motif
    `mtype`: motif_type enum
        the enum type of the given motif
    """

    __slots__ = ["residues", "mtype"]

    def __init__(self, residues, mtype):
        self.residues, self.mtype = residues, mtype


class Basepair(object):
    """
    Container for basepairs from ref_frames and dssr files produced
    by x3dna. Should not be instantiated outside x3dna.X3dna class.

    :attributes:

    `res1` : Residue
        first residue in basepair
    `res2` : Residue
        second residue in basepair
    `r` : np.array
        rotation matrix describing principle axies of orientation of basepair
    `d` : np.array
        center point of basepair
    `bp_type`: str
        x3dna dssr basepair type, see x3dna for more details

    """

    __slots__ = ["res1", "res2", "r", "d", "bp_type"]

    def __init__(self, res1, res2, r, d):
        self.res1, self.res2, self.r, self.d = res1, res2, r, d
        self.bp_type = X3dnaBPType.cDDD

    def __repr__(self):
        res1 = self.res1.chain_id + str(self.res1.num) + self.res1.i_code
        res2 = self.res2.chain_id + str(self.res2.num) + self.res2.i_code

        return "<X3DNA Basepair(%s %s %s)>" % \
            (res1, res2, self.bp_type)


class Residue(object):
    """
    Container for residues from ref_frames and dssr files produced
    by x3dna. Should not be instantiated outside x3dna.X3dna class.

    :attributes:

    `num` : int
        residue number
    `chain_id` : str
        residue chain id
    `i_code` : str
        residue insertion code

    """

    __slots__ = ["num", "chain_id", "i_code"]

    def __init__(self, num, chain_id, i_code=" "):
        self.num, self.chain_id, self.i_code = num, chain_id, i_code

    def __repr__(self):
        return "<X3DNA Residue(%s%s %s)>" % \
            (self.num, self.i_code, self.chain_id)

    def __eq__(self, res):
        if res == None:
            return 0

        return self.num == res.num and self.chain_id == res.chain_id and \
               self.i_code == res.i_code

    def __ne__(self, res):
         return not (self.num == res.num and self.chain_id == res.chain_id \
                     and self.i_code == res.i_code)

