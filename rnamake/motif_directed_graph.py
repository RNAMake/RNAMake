import motif, motif_type
import motif_type_directed_graph
from . import exceptions

class MotifDirectedGraph(motif_type_directed_graph.MotifTypeDirectedGraph):
    def __init__(self, rm):
        super(self.__class__, self).__init__()
        self._rm = rm
        self._motif_class_type = motif.Motif
        self._motif_aligner = motif.MotifAligner()

    @classmethod
    def from_str(cls, s):
        pass
        #motif_type_directed_graph.MotifTypeDirectedGraph.from_str()

    def add_motif_by_name(self, m_name, m_end_name=None, parent_index=-1,
                          parent_end_index=-1, parent_end_name=None, index=-1):

        m = self._get_motif_from_manager(m_name, m_end_name)
        return self.add_motif(m, parent_index, parent_end_index,
                              parent_end_name, index)

    def _get_motif_from_manager(self, m_name, m_end_name):
        try:
            if m_end_name is not None:
                m = self._rm.get_motif(name=m_name, end_name=m_end_name)
            else:
                m = self._rm.get_motif(name=m_name)
        except exceptions.ResourceManagerException as e:
            raise ValueError(
                "cannot add motif to graph, motif cannot be found in resource "
                "manager")
        return m

    def _split_helices_into_bp_steps(self, i, new_mdg, index_hash):
        # how many bp step motifs do we need
        name_spl = self._dg.get_node(i).name.split(".")
        if len(name_spl) == 3:
            bp_step_count = int(name_spl[2])
        else:
            bp_step_count = 1
        edges = self._dg.get_edges(i)
        h = self._rm.get_motif(name="HELIX.IDEAL")
        # add first HELIX.IDEAL motif must be connected to previous parent
        if self._dg.has_parent(i):
            parent_index = index_hash[self.get_parent_index(i)]
            new_pos = new_mdg.add_motif(h, parent_index,
                                        self.get_parent_end_index(i))
        else:
            # make sure aligned to original base pair
            self._motif_aligner.align(self.get_motif(i).get_end(0), h)
            new_pos = new_mdg.add_motif(h)

        for i in range(0, bp_step_count):
            h_new = self._motif_class_type.copy(h, new_uuid=1)
            new_pos = new_mdg.add_motif(h_new, new_pos, 1)

        return new_pos


    def get_graph_wo_ideal_helices(self):
        index_hash = {}
        new_mdg = self.__class__(self._rm)

        for i in self._dg:
            m = self._dg.get_node(i)

            # split up HELIX.IDEAL.XXX motifs into base pair steps
            if m.mtype == motif_type.HELIX and m.num_res() != 4:
                new_pos = self._split_helices_into_bp_steps(i, new_mdg, index_hash)
            else:
                new_m = self._motif_class_type.copy(m)
                if self.has_parent(i):
                    parent_index = index_hash[ self.get_parent_index(i) ]
                    new_pos = new_mdg.add_motif(new_m, parent_index,
                                                self.get_parent_end_index(i))
                else:
                    new_pos = new_mdg.add_motif(new_m)
            index_hash[i] = new_pos
        return new_mdg


    def get_secondary_structure_graph(self):
        pass

    def get_motif_state_graph(self):
        pass


    def nodes_to_pdbs(self, name="node"):
        for ni in self._dg:
          self._dg.get_node(ni).to_pdb(name+"."+str(ni)+".pdb")

