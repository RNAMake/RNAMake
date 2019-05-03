import unittest
import os
import pandas as pd
import subprocess

from rnamake.wrappers.ddg_tecto_wrapper import DDGTectoWrapper

EXE_PATH = os.environ["RNAMAKE"] + "/rnamake/lib/RNAMake/cmake/build/ddg_tecto"

class LogStatus:
    INFO = 0
    ERROR = 1

def parse_output_line(line):
    spl = line.split()
    if   spl[1] == "INFO":
        return LogStatus.INFO
    elif spl[1] == "ERROR":
        return LogStatus.ERROR

def get_cmd_output(cmd):
    s = subprocess.check_output(cmd, shell=True)
    lines = s.split("\n")
    return lines[:-1]


class DDGTectoUnittest(unittest.TestCase):

    def test_run(self):
        w = DDGTectoWrapper()
        run_opts = {
            'fseq' : "CTAGGAATCTGGAAGTACCGAGGAAACTCGGTACTTCCTGTGTCCTAG",
            'fss'  : "((((((....((((((((((((....))))))))))))....))))))",
            'cseq' : "CTAGGATATGGAAGATCCTCGGGAACGAGGATCTTCCTAAGTCCTAG",
            'css'  : "(((((((..((((((((((((....))))))))))))...)))))))"
        }
        w.run(**run_opts)
        print w.get_output()

    #def test_invalid_bp_name(self):
    #    cmd = EXE_PATH + " -pdb ttr.pdb -start_bp \"A149-A154\" -end_bp \"A222-A25\""
    #    lines = get_cmd_output(cmd)
    #    self.failUnless(parse_output_line(lines[-1]) == LogStatus.ERROR)







    """def tests(self):
        path = RESOURCE_PATH + "design_rna_tests.txt"
        df = pd.read_csv(path)

        for i, r in df.iterrows():

            cmd_opts = {
                'pdb'      : RESOURCE_PATH + r.pdb,
                'start_bp' : r.start_bp,
                'end_bp'   : r.end_bp
            }

            print r.pdb

            drw = design_rna_wrapper.DesignRNAWrapper()
            drw.run(**cmd_opts)

            self.failUnless(_out_files_exist(), "out files were not produced")
            mg = _get_motif_graph()

            self.failUnless(len(mg) > 1)
            _remove_out_files()"""



def main():
    unittest.main()

if __name__ == '__main__':
    main()

