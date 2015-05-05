import motif_library
import motif_tree
import motif
import motif_type
import settings
import sqlite3
import os
import random

class MotifLibrarySqlite(object):
    def __init__(self, libtype=None, libpath=None):
        self.motif_paths = {}
        self.connection = self._setup_sqlite3_connection(libtype, libpath)
        self.mdict = {}
        self.mtype = libtype

    def get_motif(self, name):
        if name in self.mdict:
            return self.mdict[name].copy()

        data = self.connection.execute('SELECT * from motifs WHERE name=:Name LIMIT 1',
                                       {"Name":name}).fetchone()

        if data is None:
            raise ValueError("unknown motif: " + name + "cannot get it")

        m = motif.str_to_motif(data[0])
        m.mtype = self.mtype
        self.mdict[name] = m
        return m.copy()

    def load_all(self, limit=9999):
        i = 0
        for name in self.motif_paths.iterkeys():
            i += 1
            self.get_motif(name)
            if i > limit:
                break

    def _setup_sqlite3_connection(self, libtype, libpath):
        if libtype is not None:
            try:
                name = motif_library.lib_paths[libtype]
            except:
                raise ValueError("not a valid type")
            path = settings.RESOURCES_PATH +"/motif_libraries/"+name+".db"
        else:
            path = libpath

        connection = sqlite3.connect(path)
        rows = connection.execute('SELECT name from motifs').fetchall()
        for r in rows:
            self.motif_paths[r[0]] = r[0]

        return connection

    def random_motif(self):
        return random.choice(self.motifs())

    def motifs(self):
        return self.mdict.values()

    def __contains__(self, mname):
        return mname in self.motif_paths


def build_sqlite_libraries():
    for lib_type, name in motif_library.lib_paths.iteritems():
        path = settings.RESOURCES_PATH +"/motif_libraries/"+name+".db"
        if os.path.isfile(path):
            os.remove(path)
        connection = sqlite3.connect(path)
        connection.execute("CREATE TABLE motifs(str TEXT, name TEXT, PRIMARY KEY (name))")
        try:
            mlib = motif_library.MotifLibrary(lib_type)
        except:
            continue
        mlib.load_all()
        data = []
        for m in mlib.motifs():
            s = m.to_str()
            data.append([s,m.name])
        connection.executemany("INSERT INTO motifs(str,name) VALUES (?,?) ", data)
        connection.commit()

def build_sqlite_libraries_2():
    mt = motif_tree.MotifTree()
    h_lib = motif_library.MotifLibrary(motif_type.HELIX)
    mlib = motif_library.unique_twoway_lib()
    mt.add_motif(h_lib.get_motif("HELIX.IDEAL.12"), end_index =1, end_flip=0)
    mt.level += 1

    aligned_motifs = []
    for i, m in enumerate(mlib.motifs()):
        for ei in (0, 1):
            node1 = mt.add_motif(m, end_index = ei)
            if node1 is None:
                continue

            node = mt.add_motif(h_lib.get_motif("HELIX.IDEAL.12"), end_flip=0, end_index=1)
            if node is not None:
                if m.name + "-" + str(ei) == "TWOWAY.1S72.30-1":
                    print "made it"
                    mt.write_pdbs()
                aligned_motifs.append(mt.nodes[2].motif)
                aligned_motifs[-1].name = m.name + "-" + str(ei)
                aligned_motifs[-1].ends[0].flipped=0
                aligned_motifs[-1].ends[1].flipped=0
                aligned_motifs[-1].end_to_add = ei
                aligned_motifs[-1].mtype = motif_type.TWOWAY
                mt.remove_node_level()
                continue

            other_ei = 0
            if ei == 0:
                other_ei = 1
            mt.nodes[2].motif.ends[other_ei].flip()
            node = mt.add_motif(h_lib.get_motif("HELIX.IDEAL.12"), end_flip=0, end_index=1)

            if node is not None:
                aligned_motifs.append(mt.nodes[2].motif)
                aligned_motifs[-1].ends[0].flipped=0
                aligned_motifs[-1].ends[1].flipped=0
                aligned_motifs[-1].name = m.name + "-" + str(ei)
                aligned_motifs[-1].end_to_add = ei
                aligned_motifs[-1].mtype = motif_type.TWOWAY

            mt.remove_node_level()

    path = settings.RESOURCES_PATH +"/motif_libraries/twoway_aligned.db"
    print len(aligned_motifs)
    if os.path.isfile(path):
        os.remove(path)
    connection = sqlite3.connect(path)
    connection.execute("CREATE TABLE motifs(str TEXT, name TEXT, PRIMARY KEY (name))")
    data = []
    for m in aligned_motifs:
        s = m.to_str()
        data.append([s,m.name])
    connection.executemany("INSERT INTO motifs(str,name) VALUES (?,?) ", data)
    connection.commit()

def build_sqlite_libraries_3():
    mt = motif_tree.MotifTree()
    h_lib = motif_library.MotifLibrary(motif_type.HELIX)
    mlib = motif_library.unique_twoway_lib()
    m = mlib.get_motif("TWOWAY.1S72.30")
    print m.name

    for i in range(1,13):
        mt.add_motif(h_lib.get_motif("HELIX.IDEAL."+str(i)), end_flip=0, end_index=1)
        node = mt.add_motif(m, end_index=0)
        if node:
            print i, node, mt.nodes[2].flip
        else:
            print i, node


if __name__ == '__main__':
    build_sqlite_libraries_2()
