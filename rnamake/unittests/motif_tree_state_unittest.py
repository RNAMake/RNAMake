import unittest
import rnamake.motif_tree_state
import rnamake.motif_type
import rnamake.settings
import rnamake.cluster
import rnamake.basepair
import rnamake.basic_io
import rnamake.resource_manager
import rnamake.motif
import rnamake.util
import util
import numerical
import sys
import random
import traceback

class MotifTreeStateUnittest(unittest.TestCase):

    def test_creation(self):
        #if util.UnittestState == util.UnittestType.BASIC:
        #    self.skipTest("test_creation is not a basic test")

        path = rnamake.settings.UNITTEST_PATH + "/resources/test.new.me"
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(libpath=path)
        if len(mts_lib.motif_tree_states) != 510:
            self.fail("did not load all the mts")

        clusters = rnamake.cluster.cluster_mts(mts_lib.motif_tree_states,
                                               max_distance=0.1)
        if len(clusters) != 510:
            self.fail("no all mts are unique")

    def test_parse_db_name(self):
        name = "HELIX.LE.16-0-0-0-0-1-1"
        name_elements = rnamake.motif_tree_state.parse_db_name(name)
        if name_elements.motif_name != "HELIX.LE.16":
            self.fail()

    def test_node_creation(self):
        path = rnamake.settings.UNITTEST_PATH + "/resources/test.new.me"
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(libpath=path)
        mts = mts_lib.motif_tree_states[0]
        node = rnamake.motif_tree_state.MotifTreeStateNode(mts, 0, None, 0, [0])

    def test_node_copy(self):
        path = rnamake.settings.UNITTEST_PATH + "/resources/test.new.me"
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(libpath=path)
        mts = mts_lib.motif_tree_states[0]
        node = rnamake.motif_tree_state.MotifTreeStateNode(mts, 0, None, 0, [0])
        cnode = node.copy()
        cnode.states[0].d += [0,0,-5]
        if numerical.are_points_equal(cnode.states[0].d, node.states[0].d):
            self.fail("copy was not sucessful")

    def test_aligner(self):
        return
        path = rnamake.settings.UNITTEST_PATH + \
               "/resources/test_node_align.dat"
        f = open(path)
        lines = f.readlines()
        f.close()

        path = rnamake.settings.UNITTEST_PATH + "/resources/test.new.me"
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(libpath=path)
        mts = mts_lib.motif_tree_states[0]
        node1 = rnamake.motif_tree_state.MotifTreeStateNode(mts, 0, None, 0, [0])
        node2 = rnamake.motif_tree_state.MotifTreeStateNode(mts, 0, None, 0, [0])

        aligner = rnamake.motif_tree_state.MotifTreeStateNodeAligner()
        for l in lines:
            spl = l.split(",")
            bp1 = rnamake.basepair.str_to_basepairstate(spl[0])
            bp2 = rnamake.basepair.str_to_basepairstate(spl[1])
            bp3 = rnamake.basepair.str_to_basepairstate(spl[2])
            node1.state = bp1
            node2.state = bp2
            aligner.transform_state(node1, node2)
            print node2.state.r
            print bp3.r
            exit ()

    def test_tree(self):
        mtst = rnamake.motif_tree_state.MotifTreeStateTree()

        path = rnamake.settings.UNITTEST_PATH + "/resources/test.new.me"
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(libpath=path)
        mts = mts_lib.motif_tree_states[0]
        mtst = rnamake.motif_tree_state.MotifTreeStateTree(mts)

    def test_tree_add(self):
        mtst = rnamake.motif_tree_state.MotifTreeStateTree()
        path = rnamake.settings.UNITTEST_PATH + "/resources/test.new.me"
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(libpath=path)
        mts = mts_lib.motif_tree_states[0]

        node = mtst.add_state(mts)
        if node is None or len(mtst.nodes) != 2:
            self.fail("did not add node properly")

    def test_tree_to_motiftree(self):
        mtst = rnamake.motif_tree_state.MotifTreeStateTree()
        mt = mtst.to_motiftree()
        if mt.nodes[0].motif.name != "start":
            self.fail("did not convert properly")


        path = rnamake.settings.UNITTEST_PATH + "/resources/test.new.me"
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(libpath=path)
        mts = mts_lib.motif_tree_states[0]
        mtst = rnamake.motif_tree_state.MotifTreeStateTree(mts)
        mt = mtst.to_motiftree()
        for mts in mts_lib.motif_tree_states:
            node = mtst.add_state(mts)
            if node is None:
                continue
            mt = mtst.to_motiftree()
            mt.write_pdbs()
            return

    def test_compare_mtst_to_motif_tree(self):
        if util.UnittestState == util.UnittestType.BASIC:
            self.skipTest("test_compare_mtst_to_motif_tree is not a basic test")

        mtype = rnamake.motif_type.TWOWAY
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(mtype)
        mt = None
        for j in range(1000):
            mtst = rnamake.motif_tree_state.MotifTreeStateTree()
            for i in range(100):
                mts = random.choice(mts_lib.motif_tree_states)
                mtst.add_state(mts)
            print j, len(mtst.nodes)
            try:
                mt = mtst.to_motiftree()
                mt.to_pose()
            except:
                print traceback.format_exc()
                f = open("mt.out", "w")
                f.write(mtst.to_str())
                f.close()
                exit()

    def test_compare_mtst_to_mt_features(self):
        mtype = rnamake.motif_type.TWOWAY
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(mtype)
        mtst = rnamake.motif_tree_state.MotifTreeStateTree()
        for i in range(100):
            mts = random.choice(mts_lib.motif_tree_states)
            mtst.add_state(mts)
        mt = mtst.to_motiftree()
        mt.to_pdb()
        mt_end = mt.last_node.available_ends()[0]
        # TODO figure out why numbers are slightly off

    def test_mts_to_str(self):
        path = rnamake.settings.UNITTEST_PATH + "/resources/test.new.me"
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(libpath=path)
        mts = mts_lib.motif_tree_states[0]
        s = mts.to_str()
        mts2 = rnamake.motif_tree_state.str_to_motif_tree_state(s)

    def test_tree_to_str(self):
        path = rnamake.settings.UNITTEST_PATH + "/resources/test.new.me"
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(libpath=path)
        mts = mts_lib.motif_tree_states[0]
        mtst = rnamake.motif_tree_state.MotifTreeStateTree()
        mtst.add_state(mts)
        s = mtst.to_str()
        mtst2 = rnamake.motif_tree_state.str_to_motif_tree_state_tree(s)

    def test_tree_to_str_2(self):
        return
        mtype = rnamake.motif_type.TWOWAY
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(mtype)
        mt = None
        mtst = rnamake.motif_tree_state.MotifTreeStateTree()
        for i in range(100):
            mts = random.choice(mts_lib.motif_tree_states)
            mtst.add_state(mts)

        s = mtst.to_str()
        mtst2 = rnamake.motif_tree_state.str_to_motif_tree_state_tree(s)

    def test_motif_to_state(self):
        rm = rnamake.resource_manager.ResourceManager()
        m = rm.get_motif("HELIX.IDEAL.2")
        mts = rnamake.motif_tree_state.motif_to_state(m)
        m2 = rnamake.motif.str_to_motif(mts.build_string)

    def test_compare_output(self):
        f = open("mt.out")
        l = f.readline()
        f.close()

        mtst = rnamake.motif_tree_state.str_to_motif_tree_state_tree(l)
        for i, n in enumerate(mtst.nodes):
            rnamake.basic_io.points_to_pdb("beads."+str(i)+".pdb", n.beads)

        mt = mtst.to_motiftree(sterics=1)
        for i, n in enumerate(mt.nodes):
            beads = []
            for b in n.motif.beads:
                if b.btype != 0:
                    beads.append(b.center)
            rnamake.basic_io.points_to_pdb("mt_beads."+str(i)+".pdb",beads)
        mt.write_pdbs()

    def test_bead_placement(self):
        mtype = rnamake.motif_type.TWOWAY
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(mtype)
        mt = None
        mtst = rnamake.motif_tree_state.MotifTreeStateTree()
        for mts1 in mts_lib.motif_tree_states:
            #print mts1.name
            mtst.add_state(mts1)
            for mts2 in mts_lib.motif_tree_states:
                node = mtst.add_state(mts2)
                if node is None:
                    continue
                mt = mtst.to_motiftree()
                beads = []
                for n in mtst.nodes:
                    beads.extend(n.beads)
                mt_beads = []
                for n in mt.nodes:
                    for b in n.motif.beads:
                        if b.btype == 0:
                            continue
                        mt_beads.append(b.center)
                for i in range(len(mt_beads)):
                    dist = rnamake.util.distance(mt_beads[i], beads[i])
                    if dist > 0.1:
                        print mtst.nodes[1].mts.name, mtst.nodes[2].mts.name
                        break
                mtst.remove_node(mtst.last_node)
            mtst.remove_node(mtst.last_node)

    def test_replace_state(self):
        mtype = rnamake.motif_type.TWOWAY
        mts_lib = rnamake.motif_tree_state.MotifTreeStateLibrary(mtype)
        mt = None
        mtst = rnamake.motif_tree_state.MotifTreeStateTree()
        for i in range(100):
            mts = random.choice(mts_lib.motif_tree_states)
            mtst.add_state(mts)
        mtst.to_pdb("test.pdb")
        while 1:
            mts = random.choice(mts_lib.motif_tree_states)
            node = random.choice(mtst.nodes)
            if node == mtst.nodes[0]:
                continue
            if mtst.replace_state(node, mts) == 1:
                mtst.to_pdb("test2.pdb")
                break








def main():
    unittest.main()

if __name__ == '__main__':
    # TODO make this cleaner
    if len(sys.argv) > 1:
        try:
            util.UnittestState = int(sys.argv[1])
            sys.argv.pop()
        except:
            pass
    main()
