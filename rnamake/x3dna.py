import subprocess
import settings
import os

X3DNA_BIN_PATH = settings.RESOURCES_PATH + "x3dna/bin/"
os.environ['X3DNA'] =  settings.RESOURCES_PATH + "x3dna"

# TODO figure out what operating system is being used
class X3dna(object):
    def __init__(self):
        pass

    def generate_ref_frame(self, pdb_name):
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



