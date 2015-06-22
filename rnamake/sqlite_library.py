import sqlite3
import random
import os
import settings
import motif
import motif_ensemble

class SqliteLibrary(object):
    def __init__(self):
        self.data = {}
        self.data_path = {}
        self.libnames = {}

    def _setup(self, path):
        self.connection = sqlite3.connect(path)
        rows = self.connection.execute('SELECT id from data_table').fetchall()
        for r in rows:
            self.data_path[r[0]] = r[0]

    def _get_path(self, libname):
        try:
            libpath = settings.RESOURCES_PATH  + self.libnames[libname]
        except:
            raise IOError("cannot find sqlite library")

        return libpath

    def _generate_data(self, s):
        return s

    def get(self, id):
        if id in self.data:
            return self.data[id].copy()

        if id not in self.data_path:
            raise ValueError("cannot find id in table: " + id)

        data = self.connection.execute('SELECT * from data_table WHERE id=:Id LIMIT 1',
                                       {"Id":id}).fetchone()

        if data is None:
            raise ValueError("cannot find id in table, should exist " + id)

        self.data[id] = self._generate_data(data[0])
        return self.data[id].copy()

    def get_random(self):
        return self.get(random.choice(self.data_path.keys()))

    def load_all(self, limit=99999):
        i = 0
        for id in self.data_path.keys():
            i += 1
            self.get(id)
            if i > limit:
                break

    def all(self):
        return self.data.values()

    def contains(self, id):
        if id in self.data_path:
            return 1
        else:
            return 0

class MotifSqliteLibrary(SqliteLibrary):
    def __init__(self, libname):
        super(MotifSqliteLibrary, self).__init__()
        self.libnames = self.get_libnames()
        path = self._get_path(libname)
        self._setup(path)

    def _generate_data(self, s):
        return motif.str_to_motif(s)

    @staticmethod
    def get_libnames():
        libnames = {
            "ideal_helices" :  "/motif_libraries_new/ideal_helices.db",
            "twoways"       :  "/motif_libraries_new/twoways.db"
        }

        return libnames

class MotifEnsembleSqliteLibrary(SqliteLibrary):
    def __init__(self, libname):
        super(MotifEnsembleSqliteLibrary, self).__init__()
        self.libnames = self._get_libnames()
        path = self._get_path(libname)
        self._setup(path)

    def _generate_data(self, s):
        return motif_ensemble.str_to_motif_ensemble(s)

    @staticmethod
    def get_libnames():
        libnames = {
            "bp_steps" :  "/motif_ensemble_libraries/bp_steps.db",
        }

        return libnames

class MotifStateSqliteLibrary(SqliteLibrary):

    def __init__(self, libname):
        super(MotifStateSqliteLibrary, self).__init__()
        self.libnames = self.get_libnames()
        path = self._get_path(libname)
        self._setup(path)

    @staticmethod
    def get_libnames():
        libnames = {
            "ideal_helices" :  "/motif_state_libraries/ideal_helices.db",
            "twoways"       :  "/motif_state_libraries/twoways.db"
        }

        return libnames

    def _generate_data(self, s):
        return motif.str_to_motif_state(s)


def build_sqlite_library(path, data_values, ids):
    if os.path.isfile(path):
        os.remove(path)
    connection = sqlite3.connect(path)
    connection.execute("CREATE TABLE data_table(data TEXT, id TEXT, PRIMARY KEY (id))")
    data = []
    for i, d in enumerate(data_values):
        data.append([d.to_str(),ids[i]])
    connection.executemany("INSERT INTO data_table(data,id) VALUES (?,?) ", data)
    connection.commit()