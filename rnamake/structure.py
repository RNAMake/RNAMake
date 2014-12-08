import basic_io
import pdb_parser

class Structure(object):
    def __init__(self,pdb=None):
        if pdb:
            residues = pdb_parser.parse(pdb)



