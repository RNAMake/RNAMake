import unittest
import numpy as np

from rnamake import resource_manager as rm
from rnamake import motif_factory, util, transform, motif_graph, atom, motif_tree
import is_equal


def replace_missing_phosphate_backbone(r, r_template):
    ref_frame_1 = get_res_ref_frame(r)
    ref_frame_2 = get_res_ref_frame(r_template)

    rot = util.unitarize(ref_frame_1.T.dot(ref_frame_2))
    t = transform.Transform(rot, -r_template.center())
    r_template.transform(t)
    r_template.move(r.get_atom("C4'").coords - r_template.get_atom("C4'").coords)

    atom_names = ["C5'", "O5'", "P", "OP1", "OP2"]
    new_atoms = []
    for name in atom_names:
        a = r_template.get_atom(name)
        new_atoms.append(atom.Atom(a.name, np.copy(a.coords)))

    for na in new_atoms:
        pos = r.rtype.atom_map[na.name]
        r.atoms[pos] = na


def get_res_ref_frame(r):
    beads =  r.get_beads()
    c = r.center()
    if r.name == "A" or r.name == "G":
        vec1 = util.normalize(np.array(r.get_atom("N9").coords - r.get_atom("C1'").coords))
        vec2 = util.normalize(np.array(r.get_atom("N9").coords - beads[-1].center))
        cross = np.cross(vec1, vec2)
        rot = [vec1, vec2, cross]
    else:
        vec1 = util.normalize(np.array(r.get_atom("N1").coords - r.get_atom("C1'").coords))
        vec2 = util.normalize(np.array(r.get_atom("N1").coords - beads[-1].center))
        cross = np.cross(vec1, vec2)
        rot = [vec1, vec2, cross]

    rot = np.array(rot)
    return util.unitarize(rot)


class ChainClosureUnittest(unittest.TestCase):

    def setUp(self):
        pass

    def _test_find_missing_atoms(self):
        bp_lib = rm.manager.mlibs["bp_steps"]
        bp_lib.load_all(100)
        for m in bp_lib.all():
            r = m.residues()[0]
            if r.get_atom("P") is not None:
                continue
            print m.name
            print r.name

    def _test_fix_missing_atoms(self):
        start_m = motif_factory.factory.base_motif
        m = rm.manager.get_motif(name="BP.0.27")

        r1_mis = m.residues()[0]
        r1_cor = start_m.residues()[0]

        ref_frame_1 = get_res_ref_frame(r1_mis)
        ref_frame_2 = get_res_ref_frame(r1_cor)

        r = util.unitarize(ref_frame_1.T.dot(ref_frame_2))
        t = transform.Transform(r, -r1_cor.center())
        r1_cor.transform(t)
        r1_cor.move(r1_mis.get_atom("C4'").coords - r1_cor.get_atom("C4'").coords)

        r1_mis.to_pdb("r1.pdb")
        r1_cor.to_pdb("r2.pdb")

    def test_full_fix(self):
        start_m = motif_factory.factory.base_motif
        r_template = start_m.residues()[0]

        m = rm.manager.get_motif(name="BP.0.27")
        mg = motif_graph.MotifGraph()
        """for r in m.residues():
            if r.get_atom("P") is None:
                replace_missing_phosphate_backbone(r, r_template)"""

        for i in range(10):
            mg.add_motif(m)

        #mg.to_pdb("fixed.pdb", renumber=1, close_chain=1)


    def test_tecto_fix(self):
        #m = rm.manager.get_structure("/Users/jyesselm/Documents/tecto_for_nmr/tightest.0.pdb",
        #                             "tightest")

        f = open("/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/cmake/build/bound.0.mt")
        lines = f.readlines()
        f.close()

        mt = motif_tree.motif_tree_from_topology_str(lines[0])
        mt.add_connection(1, mt.last_node().index)
        mt.to_pdb("test.pdb", renumber=1)
        mt.to_pdb("fix.pdb", renumber=1, close_chain=1)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
