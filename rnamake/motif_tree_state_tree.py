import settings
import motif_tree_state
import resource_manager
import base
import option
import util
import motif_tree
import motif

class MotifTreeStateTree(base.Base):
    def __init__(self, head_state=None, mt=None, **options):
        self.setup_options_and_constraints()
        self.options.dict_set(options)
        self.aligner = motif_tree_state.MotifTreeStateNodeAligner()
        self.clash_radius = settings.CLASH_RADIUS


        if head_state is None:
            head = self._get_default_head()
        else:
            head = self._get_head_node(head_state)

        self.nodes = [ head ]
        self.last_node = head

        if mt is not None:
            self._setup_from_mt(mt)

    def _setup_from_mt(self, mt):
        #if mt.nodes[0].motif.name != "start":
        #   self.nodes[0] = resource_manager.manager.get_state(mt.nodes[0].motif.name)

        for i, n in enumerate(mt.nodes):
            if i == 0:
                continue
            parent_index = n.parent().index
            parent_end_index = n.parent_end_index()
            end_index = n.motif_end_index(n.parent())

            mts_name = n.motif.name + "-" + str(end_index)
            mts = resource_manager.manager.get_state(mts_name)

            node = self.add_state(mts, self.nodes[parent_index],
                                  parent_end_index=parent_end_index)

            if node is None:
                raise ValueError("could not translate motiftree to motif tree state tree")

    def setup_options_and_constraints(self):
        options = { 'sterics'       : 1
                  }
        self.options = option.Options(options)
        self.constraints = {}

    def add_state(self, mts, parent=None, parent_end=None, parent_end_index=None):
        if parent is None:
            parent = self.last_node
        if parent_end is not None:
            if parent_end == parent.mts.start_index:
                return None
            parent_ends = [ parent_end ]
        elif parent_end_index is not None:
            parent_ends = [ parent.states[parent_end_index] ]
        else:
            parent_ends = parent.available_ends()

        new_node = MotifTreeStateNode(mts, len(self.nodes), parent, 0, [0])
        success = 0
        for pe in parent_ends:
            self.aligner.transform_state(pe, parent, new_node)
            self.aligner.transform_beads(new_node)
            if self.option('sterics') == 1:
                if self._steric_clash(new_node):
                    continue
            parent.add_child(new_node, pe)
            success=1
            break
        if not success:
            return None
        self.nodes.append(new_node)
        self.last_node = new_node
        return new_node

    def to_motiftree(self, **options):
        for i, n in enumerate(self.nodes):
            if i == 0:
                if n.mts.name == "start":
                    mt = motif_tree.MotifTree(**options)
                else:
                    m = motif.str_to_motif(n.mts.build_string)
                    mt = motif_tree.MotifTree(m, **options)
                continue

            m =  motif.str_to_motif(n.mts.build_string)
            parent = mt.nodes [ self.nodes.index(n.parent) ]
            parent_index = n.parent_end_index()
            mt_node = mt.add_motif(m, parent=parent, end_index=n.mts.start_index,
                                   end_flip=n.mts.flip, parent_index=parent_index)
            if mt_node is None:
                print i, n.mts.name
                raise ValueError("could not successfully convert to motiftree")

        return mt

    def _get_default_head(self):
        start_mts = motif_tree_state.ref_mts()
        start_node = MotifTreeStateNode(start_mts, 0, None, 0, [0])
        return start_node

    def _steric_clash(self, new_node, index=999):
        dist = 0
        for n in self.nodes[::-1]:
            if n.index > index:
                continue
            for b1 in new_node.beads:
                for b2 in n.beads:
                    dist = util.distance(b1, b2)
                    if dist < self.clash_radius:
                        return 1
        return 0

    def _get_head_node(self, mts):
        start_node = MotifTreeStateNode(mts, 0, None, 0, [0])
        start_node.index = 0
        return start_node

    def to_str(self):
        s = ""
        for n in self.nodes:
            s += n.to_str() + "#"
        return s

    def to_pose(self):
        return self.to_motiftree(sterics=0).to_pose()

    def to_pdb(self, fname="mtst.pdb"):
        mt = self.to_motiftree(sterics=0)
        mt.to_pdb(fname)

    def nodes_to_pdbs(self, name="node"):
        mt = self.to_motiftree()
        mt.write_pdbs(name)

    def replace_state(self, node, mts, check_clash=1):
        if len(node.mts.end_states) != len(mts.end_states):
            return 0
        for i in range(len(node.mts.end_states)):
            if node.mts.end_states[i] is None and  \
               mts.end_states[i] is not None:
                return 0
            if node.mts.end_states[i] is not None and \
               mts.end_states[i] is None:
                return 0

        old_mts = node.mts
        node.replace_mts(mts)
        open_nodes = [ node ]
        seen_nodes = []
        clash = 0
        while len(open_nodes) != 0:
            current = open_nodes.pop(0)
            seen_nodes.append(current)
            parent = current.parent
            parent_end = current.parent_end()
            if parent is None:
                continue
            self.aligner.transform_state(parent_end, parent, current)
            self.aligner.transform_beads(current)
            if check_clash:
                if self._steric_clash(current, index=current.index-1):
                    clash=1
                    break
            for c in current.children:
                if c is not None:
                    open_nodes.append(c)

        if clash:
            node.replace_mts(old_mts)
            for n in seen_nodes:
                parent = n.parent
                parent_end = n.parent_end()
                if parent is None:
                    continue
                self.aligner.transform_state(parent_end, parent, n)
                self.aligner.transform_beads(n)
            return 0
        else:
            return 1

    def remove_node(self, node=None):
        if node is None:
            node = self.last_node
        parent = node.parent
        parent_index = node.parent_end_index()
        parent.children[parent_index] = None
        self.nodes.remove(node)
        self.last_node = parent

    def remove_nodes(self, index):
        for n in self.nodes[::-1]:
            if n.index > index:
                self.remove_node(n)


class MotifTreeStateNode(object):
    def __init__(self, mts, index, parent, lib_type, children_lib_types):
        if parent is None:
            self.level = 0
        else:
            self.level = parent.level+1
        self.mts, self.parent = mts, parent
        self.lib_type, self.children_lib_types = lib_type, children_lib_types
        self.beads, self.score, self.size, self.ss_score = mts.beads, 1000, 0, 0
        self.index = index
        self.states = [None for s in self.mts.end_states]
        for i, s in enumerate(self.mts.end_states):
            if s is None:
                continue
            self.states[i] = s.copy()
        self.children = [None for s in self.states ]

    def copy(self):
        c = MotifTreeStateNode(self.mts, self.index, self.parent, self.lib_type,
                               self.children_lib_types)
        c.beads = np.copy(self.beads)
        c.score, c.size, c.ss_score = self.score, self.size, c.ss_score
        c.states = [None for s in self.states]
        c.level = self.level
        for i, s in enumerate(self.states):
            if s is None:
                continue
            c.states[i] = s.copy()
        return c

    def available_ends(self):
        states = []
        for i, state in enumerate(self.states):
            if state is None:
                continue
            if self.children[i] is None:
                states.append(state)
        return states

    def parent_end_index(self):
        if self.parent is None:
            return None
        for i, n in enumerate(self.parent.children):
            if n == self:
                return i
        raise ValueError("cannot find parent_end_index")

    def parent_end(self):
        if self.parent is None:
            return None
        for i, n in enumerate(self.parent.children):
            if n == self:
                return self.parent.states[i]
        raise ValueError("cannot find parent_end_index")

    def to_str(self):
        parent_index = -1
        if self.parent is not None:
            parent_index = self.parent.index
        s = self.mts.to_str() + "!" + str(self.parent_end_index()) + "!" + \
            str(parent_index) + "!" + str(self.lib_type) + "!"
        for state in self.states:
            if state is not None:
                s += state.to_str()
            s += "E"
        s +="!" + basic_io.point_to_str(self.children_lib_types) + "!"
        s +=basic_io.points_to_str(self.beads)
        return s

    def replace_mts(self, mts):
        self.mts = mts
        self.states = [None for s in self.mts.end_states]
        for i, s in enumerate(self.mts.end_states):
            if s is None:
                continue
            self.states[i] = s.copy()

    def add_child(self, node, end):
        i = self.states.index(end)
        self.children[i] = node

    def active_states(self):
        states = []
        for s in self.states:
            if s is not None:
                states.append(s)
        return states

    def steric_clash(self,count=9999):
        current = self.parent
        c = 0
        while current is not None:
            for b1 in self.beads:
                for b2 in current.beads:
                    dist = util.distance(b1, b2)
                    if dist < settings.CLASH_RADIUS:
                        return 1
            current = current.parent
            c += 1
            if c > count:
                break
        return 0
