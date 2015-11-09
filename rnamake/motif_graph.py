import base
import option
import settings
import motif
import util
import residue
import motif_type
import graph
import resource_manager as rm
import motif_merger
import steric_lookup


class MotifGraph(base.Base):
    def __init__(self):
        self.setup_options_and_constraints()
        self.graph = graph.GraphStatic()
        self.clash_radius = settings.CLASH_RADIUS
        self.merger = motif_merger.MotifMerger()

    def __len__(self):
        return len(self.graph)

    def copy(self):
        mg = MotifGraph()
        new_graph = self.graph.copy()
        mg.graph = new_graph
        mg.merger = self.merger.copy([n.data for n in new_graph.nodes])
        return mg

    def setup_options_and_constraints(self):
        options = { 'sterics'              : 1}

        self.options = option.Options(options)
        self.constraints = {}

    def add_motif(self, m=None, parent_index=-1, parent_end_index=-1,
                  parent_end_name=None, m_name=None, m_end_name=None):
        if m is None and m_name is not None:
            if m_end_name is not None:
                m = rm.manager.get_motif(name=m_name, end_name=m_end_name)
            else:
                m = rm.manager.get_motif(name=m_name)

        parent = self.graph.last_node
        if parent_index != -1:
            parent = self.graph.get_node(parent_index)

        if parent is None:
            m_copy = m.copy()
            m_copy.get_beads(m_copy.ends)
            m_copy.new_res_uuids()
            pos = self.graph.add_data(m_copy, -1, -1, -1, len(m_copy.ends))
            self.merger.add_motif(m_copy)
            return pos

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

            pos =  self.graph.add_data(m_added, parent.index, p, 0, len(m_added.ends))
            self.merger.add_motif(m_added, m_added.ends[0],
                                  parent.data, parent.data.ends[p])
            return pos

        return -1

    def add_motif_tree(self, mt, parent_index, parent_end_name):
        parent = self.get_node(parent_index)
        bps = parent.data.get_basepair(name=parent_end_name)
        if len(bps) == 0:
            raise ValueError("cannot find parent end in add_motif_tree")
        pei = parent.data.ends.index(bps[0])

        for i, n in enumerate(mt):
            if i == 0:
                self.add_motif(n.data, parent_index, pei)
            else:
                self.add_motif(n.data)

    def add_connection(self, i, j, i_bp_name=None, j_bp_name=None):
        node_i = self.get_node(i)
        node_j = self.get_node(j)

        node_i_indexes = []
        node_j_indexes = []
        if i_bp_name is not None:
            ei = node_i.data.get_end_index(i_bp_name)
            if not node_i.available_pos(ei):
                raise ValueError("cannot connect nodes " + str(i) + " " + str(j) +
                                 "using bp: " + i_bp_name + "as its not available")
            node_i_indexes.append(ei)
        else:
            node_i_indexes = node_i.available_children_pos()

        if j_bp_name is not None:
            ei = node_j.data.get_end_index(j_bp_name)
            if not node_j.available_pos(ei):
                raise ValueError("cannot connect nodes " + str(i) + " " + str(j) +
                                 "using bp: " + j_bp_name + "as its not available")
        else:
            node_j_indexes = node_j.available_children_pos()

        if len(node_i_indexes) > 1 or len(node_j_indexes) > 1:
            raise ValueError("cannot connect nodes " + str(i) + " " + str(j) +
                             "its unclear which ends to attach")
        if len(node_i_indexes) == 0 or len(node_j_indexes) == 0:
            raise ValueError("cannot connect nodes " + str(i) + " " + str(j) +
                             "one node has no available ends")

        self.graph.connect(i, j, node_i_indexes[0], node_j_indexes[0])
        self.merger.connect_motifs(node_i.data, node_j.data,
                                   node_i.data.ends[node_i_indexes[0]],
                                   node_j.data.ends[node_j_indexes[0]])

    def remove_motif(self, pos):
        n = self.graph.get_node(pos)
        self.merger.remove_motif(n.data)
        self.graph.remove_node(pos)

    def _add_motif_to_graph(self, m, parent=None, parent_end_index=None):
        if parent == None:
            m_copy = m.copy()
            m_copy.new_res_uuids()
            m_copy.get_beads(m_copy.ends)

            pos  =  self.graph.add_data(m_copy, -1, -1, -1, len(m_copy.ends), orphan=1)
            #self.structure.add(self.graph.get_node(pos))
            self.merger.add_motif(m_copy)
            return pos

        else:
            m_added = motif.get_aligned_motif(parent.data.ends[parent_end_index],
                                              m.ends[0],
                                              m)
            m_added.new_res_uuids()
            pos = self.graph.add_data(m_added, parent.index, parent_end_index,
                                       0, len(m_added.ends))
            #self.structure.add(self.graph.get_node(pos))
            self.merger.add_motif(m_added, m_added.ends[0],
                                  parent.data, parent.data.ends[parent_end_index])
            return pos

    def replace_ideal_helices(self):
        size = len(self.graph)
        for n in self.graph.nodes:
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
            if len(name_spl) == 3:
                count = int(name_spl[2])
            else:
                count = 1
            i = n.index
            self.remove_motif(i)

            h = rm.manager.get_motif(name="HELIX.IDEAL")
            if parent is None:
                pos = self._add_motif_to_graph(h)

            else:
                pos = self._add_motif_to_graph(h, parent, parent_end_index)

            for j in range(0, count):
                pos = self._add_motif_to_graph(h, self.graph.get_node(pos), 1)

            if other:
                self.graph.connect(pos, other.index, 1, other_end_index)
                node = self.graph.get_node(pos)
                #self.structure.connect(self.graph.get_node(pos), other)
                self.merger.connect_motifs(node.data, other.data,
                                           node.data.ends[1],
                                           other.data.ends[other_end_index])

    def secondary_structure(self):
        return self.merger.secondary_structure()

    def designable_secondary_structure(self):
        ss = self.merger.secondary_structure()

        """if self._has_ideal_helices():
            raise ValueError("make sure there are no ideal helices, "
                             "call replace_ideal_helices")"""

        for n in self.graph.nodes:
            if n.data.name != "HELIX.IDEAL":
                continue
            for r in n.data.residues():
                r_ss = ss.get_residue(uuid=r.uuid)
                if r_ss is not None:
                    r_ss.name = "N"

        return ss

    def replace_helix_sequence(self, ss):

        for n in self.graph.nodes:
            if n.data.mtype != motif_type.HELIX:
                continue
            ss_m = ss.motif(n.data.id)
            spl = ss_m.end_ids[0].split("_")
            new_name = spl[0][0]+spl[2][1]+"="+spl[0][1]+spl[2][0]
            m = rm.manager.get_motif(name=new_name)

            org_res = n.data.residues()
            new_res = m.residues()

            for i in range(len(org_res)):
                new_res[i].uuid = org_res[i].uuid

            parent = None
            parent_end_index = None
            if n.connections[0] != None:
                parent = n.connections[0].partner(n.index)
                parent_end_index = n.connections[0].end_index(parent.index)
            other = None
            other_end_index = None
            if n.connections[1] != None:
                other = n.connections[1].partner(n.index)
                other_end_index = n.connections[1].end_index(other.index)

            if parent is not None:
                m_added = motif.get_aligned_motif(parent.data.ends[parent_end_index],
                                                 m.ends[0],
                                                 m)
                n.data = m_added
            else:
                n.data = m

            self.merger.update_motif(n.data)

            if other is None:
                continue
            if other.data.mtype == motif_type.HELIX:
                continue

            m_added = motif.get_aligned_motif(n.data.ends[1],
                                              other.data.ends[0],
                                              other.data)
            other.data = m_added
            self.merger.update_motif(other.data)

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

    def write_pdbs(self, name="node"):
        for n in self.graph.nodes:
            n.data.to_pdb(name+"."+str(n.index)+".pdb")

    def to_pdb(self, name="test.pdb", renumber=-1, close_chain=0):
        return self.merger.to_pdb(name, renumber=renumber, close_chain=close_chain)

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
        return leaf_nodes

    def leafs_and_ends(self):
        leaf_nodes = []
        for n in self.graph.nodes:
            f_conn = 0
            for i, c in enumerate(n.connections):
                if n.index == 0 and i == 0 and n.data.mtype != motif_type.HAIRPIN:
                    continue
                if c is not None:
                    continue
                leaf_nodes.append([n, i])
        return leaf_nodes

    def get_node(self, i):
        return self.graph.get_node(i)

    def get_beads(self, exclude_phos=1):
        pass

    def get_steric_lookup_table(self, exclude_phos=1):
        pass

    def get_end(self, pos=-1, m_name="", m_end_name=""):
        n = None
        if   pos != -1:
            n = self.graph.get_node(pos)
        elif m_name != "":
            nodes = []
            for n in self.graph.nodes:
                if n.data.name == m_name:
                    nodes.append(n)
            if len(nodes) > 0:
                raise ValueError("cannot get end, too many motifs match name given "
                                + m_name)
            n = nodes[0]

        if m_end_name != "":
            bps = n.data.get_basepair(name=m_end_name)
            if len(bps) == 0:
                raise ValueError("found motif but " + m_end_name + "is not is a "
                                                                   "basepair in it")
            end_index = n.data.get_end_index(name=m_end_name)
            end = n.data.ends[end_index]
            #check to see if the position is available
            self.graph.get_availiable_pos(n, end_index)
            return end

        else:
            avail_pos = n.available_children_pos()
            if len(avail_pos) > 1:
                raise ValueError("too many free ends to pick one in get_end")
            if len(avail_pos) == 0:
                raise ValueError("no ends available in get_end")

            end = n.data.ends[avail_pos[0]]
            return end




































