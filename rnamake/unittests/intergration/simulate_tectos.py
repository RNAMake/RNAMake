from rnamake import base, option, motif_state_ensemble_tree
from rnamake import secondary_structure_parser, motif_type, motif_tree, exceptions
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

    def _remove_Us(self, seq):
        seq_rna = ""
        for e in seq:
            if e == 'T':
                seq_rna += 'U'
            else:
                seq_rna += e
        return seq_rna

    def _get_motifs_from_seq_and_ss(self, seq, ss):
        parser = secondary_structure_parser.SecondaryStructureParser()
        mg = parser.parse_to_motif_graph(seq, ss)

        start = 0
        motifs = []
        for n in mg.graph.nodes:
            if n.data.mtype == motif_type.TWOWAY and start == 0:
                start = 1
                continue

            if n.data.mtype == motif_type.HAIRPIN:
                break

            if not start:
                continue

            if n.data.mtype == motif_type.HELIX:
                motif = rm.manager.get_motif(end_id=n.data.end_ids[0])
                motifs.append(motif)
            elif n.data.mtype == motif_type.TWOWAY:
                try:
                    motif = rm.manager.get_motif(end_id=n.data.end_ids[0])
                    motifs.append(motifs)
                except exceptions.ResourceManagerException:
                    raise ValueError(
                        "cannot find a motif that corresponds to the sequence: " +
                        motif.sequence() + " and secondary structure: " +
                        motif.dot_bracket() + " for the simulation")
            else:
                raise ValueError(
                    "motif type: " + motif_type.type_to_str(n.data.mtype) + " is not "
                    "supported in tecto simulations.")

        return motifs

    def _get_mset(self):
        fseq = self._remove_Us(self.option('fseq'))
        fss  = self.option('fss')

        cseq = self._remove_Us(self.option('cseq'))
        css  = self.option('css')

        flow_motifs = self._get_motifs_from_seq_and_ss(fseq, fss)

        chip_motifs = self._get_motifs_from_seq_and_ss(cseq, css)

        mt = motif_tree.MotifTree()
        mt.option('sterics', 0)
        m = rm.manager.get_motif(end_id="GG_LL_CC_RR")
        mt.add_motif(m)
        mt.add_motif(m_name="GGAA_tetraloop", m_end_name="A14-A15")
        mt.add_motif(flow_motifs[1], parent_end_name="A7-A22")

        for i in range(2, len(flow_motifs)):
            mt.add_motif(flow_motifs[i])

        mt.add_motif(m_name="GAAA_tetraloop", m_end_name="A149-A154")
        mt.add_motif(chip_motifs[1], parent_end_name="A222-A251")

        for i in range(2, len(chip_motifs)):
            mt.add_motif(chip_motifs[i])

        mt.write_pdbs("new")
        mset = motif_state_ensemble_tree.MotifStateEnsembleTree(mt=mt)
        return mset

    def run(self):
        mset = self._get_mset()




