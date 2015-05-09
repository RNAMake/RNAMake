import motif_tree
import resource_manager

class MotifAssembly(object):
    def __init__(self, m=None, error_on_add=1, **options):
        self.error_on_add = error_on_add
        self.mt = motif_tree.MotifTree(m, **options)

    def add_motif(self, m=None, mname=None, end=None, end_index=None,
                  parent=None, parent_index=None, parent_end=None,
                  end_name=None, parent_end_name=None):

        if m is not None and mname is not None:
            raise ValueError("must supply motif or motif name")

        if m is None:
            m = resource_manager.manager.get_motif(mname, end_index, end_name)

        mt_node = self.mt.add_motif(m, end_index, end, 0, parent, parent_index, parent_end,
                                    end_name, parent_end_name)

        if self.error_on_add and mt_node is None:
            raise ValueError("could not add motif " + m.name)

        return mt_node
    # some wrappers for motif tree
    def write_pdbs(self, name="node"):
        self.mt.write_pdbs(name)

