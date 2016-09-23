import rnamake

def test_str_to_residue():
    path = rnamake.settings.UNITTEST_PATH + "/resources/motifs/p4p6"
    m = rnamake.motif.Motif(path)
    f = open("test_str_to_residue.dat", "w")
    for r in m.residues():
        f.write(r.to_str() + "\n")
    f.close()



test_str_to_residue()
