import warnings
import os

from rnamake import pdb_parser, exceptions, chain


path = "/Users/josephyesselman/projects/REDESIGN/resources/non-redundant-rnas"


#2ZUF: random unconnected C that prody doesnt process
#3OK4: formatting is really weird rexamine
skip_dirs = "2ZUF 3OK4".split()

f = open("seqs_in_structures.dat", "w")

dirs = []
for x in os.listdir(path):
    if os.path.isdir(path + "/" + x):
        dirs.append(x)

for d in dirs:
    if d in skip_dirs:
        continue

    pdb_path = path + "/" + d + "/" + d + ".pdb"
    print pdb_path
    try:
        with warnings.catch_warnings(record=True) as w:
            new_residues = pdb_parser.parse(pdb_path)
    except exceptions.PDBParserException:
        continue

    chains = chain.connect_residues_into_chains(new_residues)
    seq = ""
    for c in chains:
        for r in c.residues:
            seq += r.name

    f.write(d +  " " + seq + "\n")

f.close()
