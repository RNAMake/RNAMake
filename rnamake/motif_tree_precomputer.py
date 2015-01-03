import option
import base
import motif_tree
import motif_library
import motif_type
import motif_scorer
import basic_io
import util
import transformations as t
import sqlite3
import argparse

class MotifTreePrecomputerSqlite3Output(object):
    def __init__(self, name):
        self.state_connection = self._get_state_connection(name)
        self.data_connection  = self._get_data_connection(name)
        self.states, self.datas, self.names = [], [], []

    def add_data(self, state, data):
        if data.name in self.mt_names:
            return

        other_ends = ",".join([bpstate.to_str() for bpstate in state.other_ends])
        processed_state = [state.uid,state.id, state.name,
                           state.start_state.to_str(), state.end_state.to_str(),
                           state.score, state.size, other_ends]
        processed_data  = [data.name, data.build_string,
                           basic_io.beads_to_str(data.beads)]

        self.states.append(processed_state)
        self.datas.append(processed_data)
        self.names.append(data.name)

        if len(self.mt_states) > 10000:
            self._commit_data()

    def finalize(self):
        self._commit_data()

    def _commit_data(self):
        self.state_connection.executemany("INSERT INTO state(uid,id,name,       \
            start_state,end_state,score,size,other_ends) VALUES (?,?,?,?,?,?,?,?)",
            self.states)

        self.data_connection.executemany("INSERT INTO data(name,mt_str,beads)   \
                                         VALUES (?,?,?) ", self.mt_datas)

        self.state_connection.commit()
        self.data_connection.commit()

        self.datas, self.states, self.names = [], [], []

    def _get_state_connection(self, name):
        state_connection = sqlite3.connect(name+"_state.db")
        try:
            state_connection.execute("CREATE TABLE state(uid TEXT, id TEXT, \
                name TEXT, start_state TEXT, end_state TEXT, score REAL,    \
                size REAL, other_ends TEXT, PRIMARY KEY (name));");
        except:
            raise ValueError("Trying to build a new sqlite state table but  \
                             the table already exists")
        return state_connection

    def _get_data_connection(self, name):
        data_connection = sqlite3.connect(name+"_data.db")

        try:
            data_connection.execute("CREATE TABLE data(name TEXT, mt_str TEXT, beads TEXT, PRIMARY KEY (name));")

        except:
            raise ValueError("Trying to build a new sqlite data table but the table already exists")

        return data_connection


class MotifTreePrecomputerTextOutput(object):
    def __init__(self,name):
        self.f = open(name + ".new.me","w")
        self.states = []

    def add_data(self,state,data):
        for s in self.states:
            dist = s.end_state.diff(state.end_state)
            if dist < 0.1:
                return
        self.states.append(state)

        centers = []
        for i, b in enumerate(data.beads):
            if b.btype == 0: continue
            centers.append(b.center)

        self.f.write(state.name + "|" + str(state.score) + "|" + str(state.size) + \
                     "|" + str(state.flip) +  "|" + data.build_string + "|" + \
                     basic_io.points_to_str(centers) + "|" + state.end_state.to_str() + "|")

        for i, end_state in enumerate(state.other_ends):
            self.f.write(end_state.to_str() + "|" + \
                         str(state.other_end_indexes[i]) + "|")

        self.f.write("\n")

    def finalize(self):
        self.states = []


class MotifTreeStateData(object):
    def __init__(self,uid,id,name,start_state,end_state,score,size,other_ends,other_end_indexes,flip):
        self.uid = uid
        self.id = id
        self.name = name
        self.start_state = start_state
        self.end_state = end_state
        self.score = score
        self.size = size
        self.other_ends = other_ends
        self.other_end_indexes = other_end_indexes
        self.flip = flip


class MotifTreeData(object):
    def __init__(self,name,build_string,beads):
        self.name = name
        self.build_string = build_string
        self.beads = beads


class MotifTreePrecomputer(base.Base):
    def __init__(self, **options):
        self.mt = motif_tree.MotifTree()
        self.mlib = motif_library.MotifLibrary(motif_type.HELIX)
        self.scorer = motif_scorer.MotifScorer()
        self.setup_options_and_constraints()
        self.options.dict_set(options)

        if self.option('data_output') == 'sqlite3':
            self.output = MotifTreePrecomputerSqlite3Output(self.option('name'))
        else:
            self.output = MotifTreePrecomputerTextOutput(self.option('name'))

    def setup_options_and_constraints(self):
        options = { 'data_output'     : 'text',
                    'max_bps_per_end' : 11,
                    'min_bps_per_end' : 0,
                    'clash_radius'    : motif_tree.MotifTree().clash_radius,
                    'flip'            : None,
                    'motif_pos'       : None,
                    'name'            : 'test' }

        self.options = option.Options(options)
        self.constraints = {}

    def precompute_library(self, mlib):
        for m in mlib.motifs():
            self.precompute_motif(m)

    def precompute_motif(self, m):
        flip_states = (0, 1)
        motif_pos = range(len(m.ends))
        if self.option('flip') is not None:
            flip_states = self.option('flip')
        if self.option('motif_pos') is not None:
            motif_pos = self.option('motif_pos')

        for end_index in motif_pos:
            for helix_end_index in (0, 1):
                #flip the direction of the end basepair rotation
                for flip in flip_states:
                    #initial helices to add before motif
                    for hcount in range(self.option('max_bps_per_end')+1):
                        self._precompute(m, end_index, helix_end_index,
                                         flip, hcount)

    def _precompute(self, m, end_index, helix_end_index, flip, hcount):
        if hcount == 0 and helix_end_index == 1:
            return
        self.mt.remove_node_level()
        if hcount > 0:
            hmotif = self._get_helix_motif(hcount)
            self.mt.add_motif(hmotif,end_index=helix_end_index)
        motif_node = self.mt.add_motif(m, end_index=end_index, end_flip=flip)
        if motif_node is None:
            return
        avail_ends = motif_node.available_ends()
        for end in avail_ends:
            self._record_state(motif_node, end)
            self._add_end_helix(motif_node, end)

    def _get_helix_motif(self, hcount):
        if hcount == 1:
            return self.mlib.get_motif("HELIX.IDEAL")
        else:
            return self.mlib.get_motif("HELIX.IDEAL."+str(hcount))

    def _add_end_helix(self, motif_node, m_end):
        self.mt.last_node = motif_node
        for i in range(1, self.option('max_bps_per_end')+1):
            hmotif = self._get_helix_motif(i)
            node = self.mt.add_motif(hmotif, parent_end=m_end)
            if node is None:
                return
            end = self.mt.last_node.available_ends()[0]
            self._record_state(motif_node, end)
            self.mt.remove_node(node)

    def _helix_count(self, name):
        name_spl = name.split(".")
        if len(name_spl) == 2:
            return 1
        return int(name_spl[2])

    def _record_state(self, motif_node, end_bp):
        start_bp = self.mt.nodes[1].connections[0].motif_end(self.mt.nodes[1])

        pose = self.mt.get_pose()
        pose.mtype = motif_node.motif.mtype
        pose_start_bp = pose.get_basepair(bp_uuid=start_bp.uuid)[0]
        pose_end_bp   = pose.get_basepair(bp_uuid=end_bp.uuid)[0]
        start_index   = pose.ends.index(pose_start_bp)
        end_index     = pose.ends.index(pose_end_bp)
        helix_direction, start_helix_count, end_helix_count = 0, 0, 0
        if len(self.mt.nodes) > 2:
            helix_direction = self.mt.nodes[1].motif.ends.index(start_bp)
            start_helix_count = self._helix_count(self.mt.nodes[1].motif.name)
        if self.mt.last_node.motif.mtype == motif_type.HELIX:
            end_helix_count =  self._helix_count(self.mt.last_node.motif.name)
        name = motif_node.motif.name + "-" + str(helix_direction) + "-" + \
               str(start_helix_count) + "-" + str(start_index) + "-" + \
               str(end_helix_count) + "-" + str(end_index) + "-" + \
               str(motif_node.flip)
        uid = self._get_id_for_state(pose_end_bp.state(), 2.0, 0.6)
        id  = self._get_id_for_state(pose_end_bp.state(), 6.0, 0.6)
        score = self.scorer.score(pose)
        size = len(pose.residues())
        other_ends, other_end_indexes = [], []
        for i, end in enumerate(pose.ends):
            if end.uuid == pose_start_bp.uuid or end.uuid == pose_end_bp.uuid:
                continue
            other_ends.append(end.state())
            other_end_indexes.append(i)
        mt_state = MotifTreeStateData(uid, id, name, pose_start_bp.state(),
                                      pose_end_bp.state(), score, size, other_ends,
                                      other_end_indexes, motif_node.flip)
        mt_data  = MotifTreeData(name, pose.to_str(),
                                 pose.get_beads([pose_start_bp]))

        self.output.add_data(mt_state, mt_data)

    def _get_id_for_state(self, state, grid_size, euler_grid_size):
        rounded = []

        for e in state.d:
            pos = round(e / grid_size)*grid_size
            rounded.append(pos)

        angles = t.euler_from_matrix(state.r,'sxyz')
        angles = list(angles)

        for j in range(len(angles)):
            if angles[j] < 0:
                angles[j] += 3.14*2

        for a in angles:
            pos = round(a / euler_grid_size)* euler_grid_size
            rounded.append(pos)

        poss = []
        for pos in rounded:
            spos = str(pos)
            if spos == "-0.0":
                spos = "0.0"
            poss.append(spos)

        return " ".join(poss) + " "


def parse_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-t', required=True)
    parser.add_argument('--options', '-o', action="store",required=False)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    mtype = motif_type.str_to_type(args.t)
    mlib = motif_library.MotifLibrary(mtype)
    mlib.load_all()
    mtp = MotifTreePrecomputer(name=args.t,max_bps_per_end=0)
    mtp.precompute_library(mlib)


