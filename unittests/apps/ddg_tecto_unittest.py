import unittest
import os
import pandas as pd
import subprocess

from rnamake.wrappers.ddg_tecto_wrapper import DDGTectoWrapper
from rnamake.wrappers import wrapper

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
            'css'  : "(((((((..((((((((((((....))))))))))))...)))))))",
            'steps' :  1
        }
        w.run(**run_opts)
        r = w.get_results()
        self.failUnless(r.runs == 1)
        self.failUnless(r.avg == 0.0)
        self.failUnless(r.stdev == 0.0)

    def test_errors(self):
        w = DDGTectoWrapper()
        with self.assertRaises(wrapper.WrapperException) as context:
            w.run()




def main():
    unittest.main()

if __name__ == '__main__':
    main()

