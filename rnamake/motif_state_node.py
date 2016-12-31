

class NodeData(object):
    def __init__(self, ref_state):
        self.ref_state = ref_state
        self.cur_state = ref_state.copy()

    def get_end_state(self, name=None, id=None):
        return self.cur_state.get_end_state(name, id)

    def get_end_index(self, name=None, id=None):
        return self.cur_state.get_end_index(name, id)

    def name(self):
        return self.cur_state.name

    def block_end_add(self):
        return self.cur_state.block_end_add

    def end_name(self, i):
        return self.cur_state.end_names[i]

    def update_state(self):
        pass

    def uuid(self):
        return self.ref_state.uuid