import sqlite3
import random
import os
import settings
import motif
import motif_ensemble
import secondary_structure_factory
import secondary_structure
import ss_tree

class SqliteLibrary(object):
    def __init__(self):
        self.data = {}
        self.data_path = {}
        self.libnames = {}

    def _setup(self, path):
        self.connection = sqlite3.connect(path)
        cols =  self.connection.execute('PRAGMA table_info(data_table)').fetchall()
        self.keys = {}
        for c in cols:
            self.keys[c[1]] = 1

        cols = self.connection.execute('SELECT count(*) from data_table').fetchone()
        self.max_size = int(cols[0])

    def _get_path(self, libname):
        self.name = libname
        try:
            libpath = settings.RESOURCES_PATH  + self.libnames[libname]
        except:
            raise IOError("cannot find sqlite library")

        return libpath

    def _generate_data(self, s):
        return s

    def _args_to_str(self, options):
        s = ""
        for k, v in options.iteritems():
            s += k + " = " + v + ","
        return s

    def _generate_query(self, options):
        cmd = "SELECT * from data_table WHERE "
        replace = "{"
        l = len(options)
        count = 0
        for k,v in options.iteritems():
            if k not in self.keys:
                raise ValueError("attempted to use " + k + "=" + v + " in getting data from "
                                 "SqliteLibrary, this key does not exist")

            cmd += k + "='" + v + "' "
            if len(options) != count+1:
                cmd += " AND "
            count += 1

        return cmd

    def get(self, **options):

        query = self._generate_query(options)
        rows = self.connection.execute(query).fetchall()
        #if len(rows) > 1:
        #    raise ValueError("query returned too many rows, if this was expected use" \
        #                     " get_multi(): " + self._args_to_str(options) )

        if len(rows) == 0:
            raise ValueError("query returned no rows: " + self._args_to_str(options) )

        id = rows[0][-1]
        if id not in self.data:
            self.data[id] = self._generate_data(rows[0][0])

        return self.data[id].copy()

    def get_multi(self, **options):
        query = self._generate_query(options)
        rows = self.connection.execute(query).fetchall()

        if len(rows) == 0:
            raise ValueError("query returned no rows: " + self._args_to_str(options) )

        datas = []
        for r in rows:
            id = r[-1]
            if id not in self.data:
                self.data[id] = self._generate_data(r[0])
            datas.append(self.data[id].copy())
        return datas

    def get_random(self):
        id = str(random.randint(1, self.max_size-1))
        return self.get(id=id)

    def load_all(self, limit=99999):
        i = 0
        for row in self.connection.execute('SELECT * from data_table').fetchall():
            i += 1
            self.get(id=row[-1])
            if i > limit:
                break

    def all(self):
        return self.data.values()

    def contains(self, **options):
        query = self._generate_query(options)
        rows = self.connection.execute(query).fetchall()
        if len(rows) == 0:
            return 0
        else:
            return 1


class MotifSqliteLibrary(SqliteLibrary):
    def __init__(self, libname):
        super(MotifSqliteLibrary, self).__init__()
        self.libnames = self.get_libnames()
        self.ss_trees = {}
        path = self._get_path(libname)
        self._setup(path)

    def _generate_data(self, s):
        return motif.str_to_motif(s)

    @staticmethod
    def get_libnames():
        libnames = {
            "ideal_helices"   : "/motif_libraries_new/ideal_helices.db",
            "ideal_helices_reversed" :  "/motif_libraries_new/ideal_helices_reversed.db",
            "twoway"          : "/motif_libraries_new/twoway.db",
            "tcontact"        : "/motif_libraries_new/tcontact.db",
            "hairpin"         : "/motif_libraries_new/hairpin.db",
            "nway"            : "/motif_libraries_new/nway.db",
            "unique_twoway"   : "/motif_libraries_new/unique_twoway.db",
            "bp_steps"        : "/motif_libraries_new/bp_steps.db",
        }

        return libnames

    def get_best_matches(self, new_id):
        #if self.contains(end_id=new_id):
        #    return self.get(end_id=new_id)

        if len(self.ss_trees) == 0:
            self.load_all()
            for m in self.all():
                sstree = secondary_structure_factory.ss_id_to_ss_tree(m.end_ids[0])
                self.ss_trees[m.end_ids[0]] = sstree

        new_ss_tree = secondary_structure_factory.ss_id_to_ss_tree(new_id)
        best_score = 10000
        best_id = None
        matches = []
        for id, sstree in self.ss_trees.iteritems():
            score = ss_tree.compare_ss_tree(new_ss_tree, sstree)
            if score < best_score:
                best_score = score
                best_id = id
            matches.append([id, score])

        if best_score == 10000:
            raise  ValueError("get_best_match failed in MotifSSIDSqliteLibrary")

        matches.sort(key=lambda x: x[1])

        motifs = []
        for i in range(0, 9):
            motifs.append(self.get(end_id=matches[i][0]))
        return motifs


class MotifSSIDSqliteLibrary(SqliteLibrary):
    def __init__(self, libname):
        super(MotifSSIDSqliteLibrary, self).__init__()
        self.libnames = self.get_libnames()
        self.libname = libname
        self.ss_trees = {}
        path = self._get_path(libname)
        self._setup(path)

    def _generate_data(self, s):
        return motif.str_to_motif(s)

    @staticmethod
    def get_libnames():
        libnames = {
            "bp_steps"     : "/motif_libraries_new/bp_steps.db",
            "twoway"       : "/motif_libraries_new/ss_twoway.db",
            "tcontact"     : "/motif_libraries_new/ss_tcontact.db",
            "hairpin"      : "/motif_libraries_new/ss_hairpin.db",
            "nway"         : "/motif_libraries_new/ss_nway.db",
        }

        return libnames






        #check end basepairs

    def get_by_topology(self, top):
        scaled_top = []
        motifs = []
        for t in top:
            scaled_top.append(t + 2)

        for id in self.data_path:
            spl = id.split("_")
            sizes = []
            for i in range(0, len(spl)-1, 2):
                sizes.append(len(spl[i+1]))

            mismatch = 0
            for i in range(len(sizes)):
                if sizes[i] != scaled_top[i]:
                    mismatch = 1
                    continue

            if not mismatch:
                motifs.append(self.get(id))

        return motifs


class MotifEnsembleSqliteLibrary(SqliteLibrary):
    def __init__(self, libname):
        super(MotifEnsembleSqliteLibrary, self).__init__()
        self.libnames = self.get_libnames()
        path = self._get_path(libname)
        self._setup(path)

    def _generate_data(self, s):
        return motif_ensemble.str_to_motif_ensemble(s)

    @staticmethod
    def get_libnames():
        libnames = {
            "bp_steps" :  "/motif_ensemble_libraries/bp_steps.db",
            "twoway"   :  "/motif_ensemble_libraries/twoway.db",
            "nway"     :  "/motif_ensemble_libraries/nway.db",
            "tcontact" :  "/motif_ensemble_libraries/tcontact.db",
            "hairpin"  :  "/motif_ensemble_libraries/hairpin.db",
            "twoway_clusters" : "motif_ensemble_libraries/twoway_clusters.db"
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
            "twoway"        :  "/motif_state_libraries/twoway.db",
            "unique_twoway" :  "/motif_state_libraries/unique_twoway.db",
            "bp_steps"   :  "/motif_state_libraries/bp_steps.db",
            "twoway"     :  "/motif_state_libraries/twoway.db",
            "nway"          :  "/motif_state_libraries/nway.db",

        }

        return libnames

    def _generate_data(self, s):
        return motif.str_to_motif_state(s)


class MotifStateEnsembleSqliteLibrary(SqliteLibrary):
    def __init__(self, libname):
        super(MotifStateEnsembleSqliteLibrary, self).__init__()
        self.libnames = self.get_libnames()
        path = self._get_path(libname)
        self._setup(path)

    def _generate_data(self, s):
        return motif_ensemble.str_to_motif_state_ensemble(s)

    @staticmethod
    def get_libnames():
        libnames = {
            "bp_steps" :  "/motif_state_ensemble_libraries/bp_steps.db",
            "twoway"   :  "/motif_state_ensemble_libraries/twoway.db",
            "nway"     :  "/motif_state_ensemble_libraries/nway.db",
            "tcontact" :  "/motif_state_ensemble_libraries/tcontact.db",
            "hairpin"  :  "/motif_state_ensemble_libraries/hairpin.db",
        }

        return libnames


class MotifClusterSqliteLibrary(SqliteLibrary):

    def __init__(self, libname):
        super(MotifClusterSqliteLibrary, self).__init__()
        self.libnames = self.get_libnames()
        path = self._get_path(libname)
        self._setup(path)

    @staticmethod
    def get_libnames():
        libnames = {
            "twoways"      :  "/motif_state_libraries/twoways_clusters.db"
        }
        return libnames

    def _generate_data(self, s):
        return motif.str_to_motif_array(s)


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


def build_sqlite_library_2(path, data, keys, primary_key):
    if os.path.isfile(path):
        os.remove(path)
    connection = sqlite3.connect(path)
    insert_str = "INSERT INTO data_table("
    table_str = "CREATE TABLE data_table("
    for key in keys:
        table_str += key + " TEXT,"
        insert_str += key
        if key != keys[-1]:
            insert_str += ","
    table_str += "PRIMARY KEY ("+primary_key+"))"
    insert_str += ") VALUES ("
    for key in keys:
        insert_str += "?"
        if key != keys[-1]:
            insert_str += ","
    insert_str += ")"

    connection.execute(table_str)
    connection.executemany(insert_str, data)
    connection.commit()














