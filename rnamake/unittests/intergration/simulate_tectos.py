from rnamake import base, option, motif_state_ensemble_tree
from rnamake import secondary_structure_parser, motif_type, motif_tree
from rnamake import resource_manager as rm

class SimulateTectos(base.Base):
    def __init__(self, **options):
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
        parser = secondary_structure_parser.SecondaryStructureParser()
        mg = parser.parse_to_motif_graph(seq, ss)

        start = 0
        names = []
        for n in mg.graph.nodes:
            if n.data.mtype == motif_type.TWOWAY:
                start = 1
                continue

            if n.data.mtype == motif_type.HAIRPIN:
                break

            if not start:
                continue

            name = self._get_bp_name_from_seq(n.data.sequence())
            new_name = ""
            for e in name:
                if e == "T":
                    new_name += "U"
                else:
                    new_name += e
            names.append(new_name)
        return names


    def _get_mset(self):
        flow_motif_names = self._get_m_names_from_seq_and_ss(self.option('fseq'),
                                                             self.option('fss'))

        chip_motif_names = self._get_m_names_from_seq_and_ss(self.option('cseq'),
                                                             self.option('css'))

        mt = motif_tree.MotifTree()
        mt.option('sterics', 0)
        mt.add_motif(m_name="GC=GC")
        mt.add_motif(m_name="GGAA_tetraloop", m_end_name="A14-A15")
        mt.add_motif(m_name=flow_motif_names[1], parent_end_name="A22-A7")

        for i in range(2, len(flow_motif_names)):
            mt.add_motif(m_name=flow_motif_names[i])

        mt.add_motif(m_name="GAAA_tetraloop", m_end_name="A149-A154")
        mt.add_motif(m_name=chip_motif_names[1], parent_end_name="A222-A251")

        for i in range(2, len(chip_motif_names)):
            mt.add_motif(m_name=chip_motif_names[i])

        mt.write_pdbs("new")
        mset = motif_state_ensemble_tree.MotifStateEnsembleTree(mt=mt)
        return mset



