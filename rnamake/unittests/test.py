import rnamake.resource_manager as rm
import rnamake.pose_factory as pf
import rnamake.motif_tree as motif_tree

rm.manager.add_motif('mfP2ab_short.pdb')
#rm.manager.add_motif('mfPK-coot-8.pdb')

mt = motif_tree.MotifTree(sterics=0)

#m = rm.manager.get_motif(name='mfP2ab_short', end_name='A1-A36')
#mt.add_motif(m)
p = pf.factory.pose_from_file('mfPK-coot-8.pdb')
mt.add_motif(p)
mt.add_motif(rm.manager.get_motif(name='mfP2ab_short', end_name='A17-A26'))
mt.write_pdbs()
mt.to_pdb("test.pdb")