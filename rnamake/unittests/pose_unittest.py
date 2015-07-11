import unittest
import rnamake.pose
import rnamake.settings as settings
import rnamake.motif_library as motif_library
import rnamake.motif_type as motif_type
import rnamake.motif_tree as motif_tree
import random
import rnamake.eternabot.sequence_designer as sequence_designer
import rnamake.pose_factory as pf


def get_unique_twoway_mlib():
    path = settings.MOTIF_DIRS + "two_ways/unique_7.dat"
    mlib = motif_library.MotifLibrary(libfile=path)
    mlib.load_all()
    return mlib

def get_ideal_helix_mlib():
    mlib = motif_library.MotifLibrary(motif_type.HELIX)
    for i in range(1,21):
        mlib.get_motif("HELIX.LE."+str(i))
    return mlib

def get_twoway_helix_motiftree(size=10):
    mlib1 = get_unique_twoway_mlib()
    mlib2 = get_ideal_helix_mlib()
    mlibs = [mlib2, mlib1]
    mt = motif_tree.MotifTree()
    count = 0
    pos = 0
    i = 0
    while i < size:
        if i % 2 == 0:
            pos = 0
        else:
            pos = 1

        m = random.choice(mlibs[pos].motifs())
        node = mt.add_motif(m)
        if node:
            i += 1
        count += 1
        if count > 1000:
            break
    return mt


class PoseUnittest(unittest.TestCase):

    def test_creation(self):
        p = pf.factory.pose_from_file("resources/motifs/p4p6")

    def test_designable_sequence(self):
        mt = get_twoway_helix_motiftree()
        p = mt.to_pose()
        seq = p.designable_sequence()
        ss = p.secondary_structure()
        #print seq
        #print ss
        designer = sequence_designer.SequenceDesigner()
        results = designer.design(ss, seq)
        #print results[0]['end'][0]

    def test_optimized_sequence(self):
        mt = get_twoway_helix_motiftree()
        p = mt.to_pose()
        seq = p.optimized_sequence()

    def test_motifs(self):
        p = pf.factory.pose_from_file("resources/motifs/p4p6")
        twoways = p.motifs(motif_type.TWOWAY)
        for i, m in enumerate(twoways):
            m.to_pdb("motif."+str(i)+".pdb")





def main():
    unittest.main()

if __name__ == '__main__':
    main()
