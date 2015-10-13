import base
import option
import settings
import motif
import util
import residue
import motif_type
import motif_tree_merger
import graph
import resource_manager as rm


class MotifGraph(base.Base):
    def __init__(self):
        self.setup_options_and_constraints()
        self.graph = graph.GraphStatic()
        self.clash_radius = settings.CLASH_RADIUS
        self.merger = motif_tree_merger.MotifTreeMerger()

    def copy(self):
        mg = MotifGraph()
        new_graph = self.graph.copy()
        mg.graph = new_graph
        return mg

    def setup_options_and_constraints(self):
        options = { 'sterics'              : 1}

        self.options = option.Options(options)
        self.constraints = {}

    def add_motif(self, m=None, parent_index=-1, parent_end_index=-1,
                  parent_end_name=None):
        parent = self.graph.last_node
        if parent_index != -1:
            parent = self.graph.get_node(parent_index)

        if parent is None:
            m_copy = m.copy()
            m_copy.get_beads(m_copy.ends)
            return self.graph.add_data(m_copy, -1, -1, -1, len(m_copy.ends))

        if parent_end_name is not None:
            parent_end = parent.data.get_basepair(name=parent_end_name)[0]
            parent_end_index = parent.data.ends.index(parent_end)

        avail_pos = self.graph.get_availiable_pos(parent, parent_end_index)

        for p in avail_pos:
            if p == parent.data.block_end_add:
                continue

            m_added = motif.get_aligned_motif(parent.data.ends[p], m.ends[0], m)
            if self.option('sterics'):
                if self._steric_clash(m_added):
                    continue

            m_added.new_res_uuids()

            return self.graph.add_data(m_added, parent.index, p, 0, len(m_added.ends))
        #self._update_beads(parent, new_node)

        return -1

    def _add_motif_to_graph(self, m, parent=None, parent_end_index=None):
        if parent == None:
            m_copy = m.copy()
            m_copy.get_beads(m_copy.ends)
            m_copy.new_res_uuids()
            return self.graph.add_data(m_copy, -1, -1, -1, len(m_copy.ends), orphan=1)

        else:
            m_added = motif.get_aligned_motif(parent.data.ends[parent_end_index],
                                              m.ends[0],
                                              m)
            m_added.new_res_uuids()
            return self.graph.add_data(m_added, parent.index, parent_end_index,
                                       0, len(m_added.ends))

    def replace_ideal_helices(self):
        size = len(self.graph)
        for i in range(0, size):
            n = self.graph.get_node(i)
            if n.data.mtype != motif_type.HELIX:
                continue
            if len(n.data.residues()) == 4:
                continue

            parent = None
            parent_end_index = None
            other = None
            other_end_index = None
            if n.connections[0] != None:
                parent = n.connections[0].partner(n.index)
                parent_end_index = n.connections[0].end_index(parent.index)
            if n.connections[1] != None:
                other = n.connections[1].partner(n.index)
                other_end_index = n.connections[1].end_index(other.index)


            name_spl = n.data.name.split(".")
            count = int(name_spl[2])
            self.graph.remove_node(i)
            h = rm.manager.get_motif(name="HELIX.IDEAL")
            if parent is None:
                pos = self._add_motif_to_graph(h)
            else:
                pos = self._add_motif_to_graph(h, parent, parent_end_index)

            for j in range(0, count):
                pos = self._add_motif_to_graph(h, self.graph.get_node(pos), 1)

            if other:
                self.graph.connect(pos, other.index, 1, other_end_index)

    def designable_secondary_structure(self, leaf=None, leaf_pos=None):
        if self._has_ideal_helices():
            raise ValueError("make sure there are no ideal helices, "
                             "call replace_ideal_helices")

        start = self.graph.oldest_node()
        if leaf is not None:
            start = leaf
        if leaf_pos is not None:
            start = self.graph.get_node(leaf_pos)

        ss = self.merger.get_secondary_structure(self.graph, start)

        for r in ss.residues():
            for n in self.graph:
                if n.data.name != "HELIX.IDEAL":
                    continue
                if n.data.get_residue(uuid=r.uuid) != None:
                    r.name = "N"
                    break

        return ss

    def replace_helix_sequence(self, ss):
        for n in self.graph:
            if n.data.mtype != motif_type.HELIX:
                continue
            for r in n.data.residues():
                try:
                    print ss.get_residue(uuid=r.uuid).name
                except:
                    pass

            exit()



    def _steric_clash(self, m):
        beads = m.beads
        for n in self.graph:
            for c1 in n.data.beads:
                for c2 in beads:
                    if c1.btype == residue.BeadType.PHOS or \
                       c2.btype == residue.BeadType.PHOS:
                        continue
                    dist = util.distance(c1.center, c2.center)
                    if dist < self.clash_radius:
                        return 1
        return 0

    def _has_ideal_helices(self):
        for n in self.graph:
            if n.data.mtype == motif_type.HELIX and len(n.data.residues()) > 4:
                return 1
        return 0

    def write_pdbs(self,name="node"):
        for n in self.graph.nodes:
            n.data.to_pdb(name+"."+str(n.index)+".pdb")

    def last_node(self):
        return self.graph.last_node

    def leafs(self):
        leaf_nodes = []
        for n in self.graph:
            f_conn = 0
            for c in n.connections:
                if c is None:
                    f_conn += 1
            if f_conn == 0:
                continue
            leaf_nodes.append(n)
        return nodes


