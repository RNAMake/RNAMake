from rnamake import base, option, motif_state_ensemble_tree, settings
from rnamake import secondary_structure_parser, motif_type, motif_tree, ss_tree
from rnamake import resource_manager as rm

class SimulateTectos(base.Base):
    def __init__(self, **options):
        base_path = settings.base_dir + "/lib/RNAMake/apps/simulate_tectos/resources/"
        rm.manager.add_motif(base_path+"GAAA_tetraloop")
        rm.manager.add_motif(base_path+"GGAA_tetraloop")

        self.setup_options_and_constraints()
        self.options.dict_set(options)
        self.mset = self._get_mset()

    def setup_options_and_constraints(self):
        options = { 'fseq'   : 'CTAGGAATCTGGAAGTACCGAGGAAACTCGGTACTTCCTGTGTCCTAG',
                    'fss'    : '((((((....((((((((((((....))))))))))))....))))))',
                    'cseq'   : 'CTAGGATATGGAAGATCCTCGGGAACGAGGATCTTCCTAAGTCCTAG',
                    'css'    : '(((((((..((((((((((((....))))))))))))...)))))))'}
        self.options = option.Options(options)

    def _get_bp_name_from_seq(self, seq):
        name = seq[0]+seq[4]+"="+seq[1]+seq[3]
        return name

    def _get_m_names_from_seq_and_ss(self, seq, ss):
        sstree = ss_tree.SS_Tree(seq, ss)

        names = []
        current = None
        for n in sstree:
            if n.data.type == ss_tree.SS_Type.SS_BULGE:
                current = n
                break

        current = current.children[0]
        required_nodes = []
        while current.data.type != ss_tree.SS_Type.SS_HAIRPIN:
            required_nodes.append(current)
            current = current.children[0]

        for i in range(1, len(required_nodes)):
            seq1 = required_nodes[i-1].data.sequence()
            seq2 = required_nodes[i].data.sequence()
            motif_name = seq1[0]+seq1[2]+"="+seq2[0]+seq2[2]
            motif_name_rna = ""
            for e in motif_name:
                if e == "T":
                    motif_name_rna += "U"
                else:
                    motif_name_rna += e
            names.append(motif_name_rna)

        return names


    def _get_mset(self):
        flow_motif_names = self._get_m_names_from_seq_and_ss(self.option('fseq'),
                                                             self.option('fss'))

        chip_motif_names = self._get_m_names_from_seq_and_ss(self.option('cseq'),
                                                             self.option('css'))

        mt = motif_tree.MotifTree()
        mt.option('sterics', 0)
        m = rm.manager.get_motif(name="GC=GC")
        mt.add_motif(m)
        m = rm.manager.get_motif(name="GGAA_tetraloop", end_name="A14-A15")
        mt.add_motif(m)
        m = rm.manager.get_motif(name=flow_motif_names[1])
        mt.add_motif(m, parent_end_name="A22-A7")

        for i in range(2, len(flow_motif_names)):
            m = rm.manager.get_motif(name=flow_motif_names[i])
            mt.add_motif(m)

        m = rm.manager.get_motif(name="GAAA_tetraloop", end_name="A149-A154")
        mt.add_motif(m)
        m = rm.manager.get_motif(name=chip_motif_names[1])
        mt.add_motif(m, parent_end_name="A222-A251")

        for i in range(2, len(chip_motif_names)):
            m = rm.manager.get_motif(name=chip_motif_names[i])
            mt.add_motif(m)

        mset = motif_state_ensemble_tree.MotifStateEnsembleTree(mt=mt)
        return mset




