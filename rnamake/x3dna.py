import subprocess
import settings
import os
import re
import numpy as np
import util

X3DNA_BIN_PATH = settings.RESOURCES_PATH + "x3dna/bin/"
os.environ['X3DNA'] =  settings.RESOURCES_PATH + "x3dna"

# TODO figure out what operating system is being used
# TODO create at enum type for basepair types instead of strings
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
        # to get the basepairing information from a pdb, for example test.pdb
        >>>x3dna = x3dna.X3dna()
        >>>basepairs = x3dna.get_basepairs("test")

    """

    def __init__(self):
        pass

    def generate_ref_frame(self, pdb_name):
        """
        runs x3dna/find_pair and x3dna/analyze to get the basepairing
        information of a given pdb and it is updated to a ref_frames.dat
        file which can be parsed by get_basepair_info if necessary

        :param pdb_name: the pdb you want to run without the .pdb extension
        :type pdb_name: str

        .. code-block:: python
            # assumes there is a file "test.pdb" in current dir
            # after running should be a "ref_frames.dat" in current dir
            >>>x3dna = rnamake.x3dna.X3dna()
            >>>x3dna.generate_ref_frame("test")

        """
        if not os.path.isfile(pdb_name + ".pdb"):
            raise IOError(pdb_name + ".pdb is not found cannot generate "+ \
                          "ref_frames.dat file")

        find_pair_path = X3DNA_BIN_PATH + "find_pair "
        analyze_path = X3DNA_BIN_PATH + "analyze "

        result = \
        subprocess.call(find_pair_path + pdb_name + ".pdb 2> /dev/null "+
                        "stdout | " + analyze_path + "stdin >& /dev/null",
                        shell=True)

        if result != 0:
            raise SystemError("find_pair did not run correctly")

        files = ("auxiliary.par,bestpairs.pdb,bp_helical.par,bp_order.dat,"+\
                "bp_step.par,cf_7methods.par,col_chains.scr,col_helices.scr,"+\
                "hel_regions.pdb,hstacking.pdb,poc_haxis.r3d,stacking.pdb").split(",")

        name_spl = pdb_name.split("/")
        files.append(name_spl[-1]+".out")

        for f in files:
            try:
                os.remove(f)
            except:
                pass

    def generate_dssr_file(self, pdb_name):
        """
        runs x3dna's dssr which assigns the secondary structure and each of
        the motifs in a RNA structure

        :param pdb_name: the pdb you want to run without the .pdb extension
        :type pdb_name: str

        """
        if not os.path.isfile(pdb_name + ".pdb"):
            raise IOError(pdb_name + ".pdb is not found cannot generate "+ \
                          "dssr file")

        name_spl = pdb_name.split("/")
        dssr_path = X3DNA_BIN_PATH + "x3dna-dssr "
        subprocess.call(dssr_path + "-i="+pdb_name+".pdb -o="+name_spl[-1]+ \
                        "_dssr.out --non-pair >& /dev/null", shell=True)
        files = ("dssr-2ndstrs.ct,dssr-2ndstrs.dbn,dssr-helices.pdb,"+ \
                "dssr-pairs.pdb,dssr-stems.pdb,hel_regions.pdb,hstacking.pdb,"+ \
                "poc_haxis.r3d,stacking.pdb,dssr-torsions.dat,dssr-Kturns.pdb,"+ \
                "dssr-multiplets.pdb,dssr-hairpins.pdb,dssr-Aminors.pdb").split(",")
        for f in files:
            try:
                os.remove(f)
            except:
                pass

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
                    raise ValueError("did not parse 3DNA ref frame correctly")
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

    def _get_ref_frame_path(self, pdb_name):
        """
        determines were the ref_frames.dat file generated by x3dna is located
        checks where the pdb is where would be the it would be if its a motif
        directory. Also checks local directory and if does not exist in either
        location it creates a new one using generate_ref_frame. This is a
        helper function do not call directly!

        :param pdb_name: the pdb you want to run without the .pdb extension
        :type pdb_name: str

        """
        base_dir = util.base_dir(pdb_name)
        ref_frames_path = None
        if   os.path.isfile(base_dir + "/ref_frames.dat"):
            ref_frames_path = base_dir + "/ref_frames.dat"
        elif os.path.isfile("ref_frames.dat"):
            ref_frames_path = "ref_frames.dat"
        else:
            self.generate_ref_frame(pdb_name)
            ref_frames_path = "ref_frames.dat"

        return ref_frames_path

    def _get_dssr_file_path(self, pdb_path):
        """
        determines were the pdb_dssr.out file generated by x3dna is located
        checks where the pdb is where would be the it would be if its a motif
        directory. Also checks local directory and if does not exist in either
        location it creates a new one using generate_dssr_file. This is a
        helper function do not call directly!

        :param pdb_name: the pdb you want to run without the .pdb extension
        :type pdb_name: str
        """
        base_dir = util.base_dir(pdb_path)
        pdb_name = util.filename(pdb_path)
        dssr_name = pdb_name + "_dssr.out"

        dssr_file_path = None
        if   os.path.isfile(base_dir + "/" + dssr_name):
            dssr_file_path = base_dir + "/" + dssr_name
        elif os.path.isfile(dssr_name):
            dssr_file_path = dssr_name
        else:
            self.generate_dssr_file(pdb_path)
            dssr_file_path = dssr_name

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
        i_code = ""

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

        :param pdb_path: path to the pdb file with out .pdb extension
        :type pdb_path: str
        """
        ref_frames_path = self._get_ref_frame_path(pdb_path)
        dssr_file_path = self._get_dssr_file_path(pdb_path)

        basepairs = self._parse_ref_frame_file(ref_frames_path)
        section = self._divide_dssr_file_into_sections(dssr_file_path)['bases']

        for l in section:
            spl = re.split("\s+",l)
            if len(spl) < 6:
                continue
            if not spl[1].isdigit():
                continue
            res1 = self._parse_dssr_res_str(spl[2])
            res2 = self._parse_dssr_res_str(spl[3])
            bp_type = spl[8]

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


class Basepair(object):

    __slots__ = ["res1", "res2", "r", "d", "bp_type"]

    def __init__(self, res1, res2, r, d):
        self.res1, self.res2, self.r, self.d = res1, res2, r, d
        self.bp_type = "c..."

class Residue(object):

    __slots__ = ["num", "chain_id", "i_code"]

    def __init__(self, num, chain_id, i_code=""):
        self.num, self.chain_id, self.i_code = num, chain_id, i_code

    def __repr__(self):
        return "<X3DNA Residue(%s%s %s)>" % \
            (self.num, self.i_code, self.chain_id)

    def __eq__(self, res):
        return self.num == res.num and self.chain_id == res.chain_id and \
               self.i_code == res.i_code

    def __ne__(self, res):
         return not (self.num == res.num and self.chain_id == res.chain_id \
                     and self.i_code == res.i_code)


