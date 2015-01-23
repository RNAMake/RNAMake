import motif_library
import motif
import settings
import sqlite3
import os

class MotifLibrarySqlite(object):
    def __init__(self, libtype):
        self.connection = self._setup_sqlite3_connection(libtype)
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
        rows = self.connection.execute('SELECT * from motifs').fetchall()
        for i, r in enumerate(rows):
            m = motif.str_to_motif(r[0])
            m.mtype = self.mtype
            self.mdict[r[1]] = m
            if i > limit:
                break

    def _setup_sqlite3_connection(self, libtype):
        try:
            name = motif_library.lib_paths[libtype]
        except:
            raise ValueError("not a valid type")
        path = settings.RESOURCES_PATH +"/motif_libraries/"+name+".db"

        connection = sqlite3.connect(path)
        return connection


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



if __name__ == '__main__':
    build_sqlite_libraries()
