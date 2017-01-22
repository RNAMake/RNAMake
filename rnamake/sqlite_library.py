import sqlite3
import random
import os

import settings
import exceptions
import motif
import motif_ensemble
import residue_type

class SqliteLibrary(object):
    """
    A wrapper for sqlite3 library, makes it simpler to do repeated sqlite3
    calls with a set database file. This is a base class and should not be
    called directly.

    :attributes:
    `data` : Dictionary of stored data
        Holds data already queryed from database file so dont have to do an
        extra query.

    """

    def __init__(self, libname):
        self.data = {}
        self._setup(libname)

    @staticmethod
    def get_libnames():
        """
        get dictionary of valid registered library names / paths. Only library
        names registered in this corresponding function in children classes
        can be loaded.

        :return: Dictionary of names / pathes for each valid sqlite3 database
        :rtype: dict
        """
        return {}

    def _setup(self, libname):
        """
        Initates sqlite3 connection with database and ensures everything is
        working properly.

        :param path: path to sqlite3 database
        :type path: string
        :return: None
        """

        path = self._get_path(libname)

        if not os.path.isfile(path):
            raise exceptions.SqliteLibraryException(
                "path: " + path + " should exist but doesn't, either it was " +
                "deleted by accident or not commited on this branch")

        self.connection = sqlite3.connect(path)
        cols =  self.connection.execute('PRAGMA table_info(data_table)').fetchall()
        self.keys = {}
        for c in cols:
            self.keys[c[1]] = 1

        cols = self.connection.execute('SELECT count(*) from data_table').fetchone()
        self.max_size = int(cols[0])

        # make sure there is something in the library
        if self.max_size == 0:
            raise exceptions.SqliteLibraryException(
                "opened database file: " + path + " but there are no rows!")

    def _get_path(self, libname):
        """
        helper method to check to make sure given name is a valid sqlite
        library.

        :param libname: name of sqlite library (twoway, nway etc)
        :type libname: str

        :return: full path of library
        :rtype: str
        """

        self.name = libname
        try:
            libpath = settings.RESOURCES_PATH  + self.get_libnames()[libname]
        except:
            raise exceptions.SqliteLibraryException(
                "libname: " + libname + " is not a valid sqlite_library. " +
                "options are: " + " ".join(self.get_libnames().keys()))

        return libpath

    def _generate_data(self, s):
        """
        do not call directly. Needs to be overloaded to un-stringify data back
        into an object.

        :param s: string version of an object stored in database
        """
        return s

    def _args_to_str(self, options):
        """
        returns string verision of options supplied used in sqlite select.
        This is mainly used to report an error if a select does not yeild
        results.

        :param options: options for SELECT query.
        :type options: dict
        :return:
        """
        s = ""
        for k, v in options.iteritems():
            s += k + " = " + v + ","
        return s

    def _generate_query(self, options):
        """
        Generates a SELECT query to retrieve rows from the sqlite3 database.

        :param options: dictionary of variables used in selection
        :type options:

        :return: sqlite3 select command
        :rtype: string
        """

        cmd = "SELECT * from data_table WHERE "
        l = len(options)
        count = 0
        for k, v in options.iteritems():
            if k not in self.keys:
                raise exceptions.SqliteLibraryException(
                    "attempted to use " + k + "=" + v + " in getting data from"
                    " SqliteLibrary, this column does not exist in database")

            cmd += k + "='" + v + "' "
            if len(options) != count+1:
                cmd += " AND "
            count += 1

        return cmd

    def get(self, **options):
        """
        Gets row for sqlite3 database with specific variables specified in
        options. Will at most return 1 item even if multiple items meet
        selection criteria.

        :param options: sqlite3 select columns and values
        :return: unstringified data from database
        """

        query = self._generate_query(options)
        rows = self.connection.execute(query).fetchall()

        if len(rows) == 0:
            raise exceptions.SqliteLibraryException(
                "query returned no rows: " + self._args_to_str(options) )

        id = rows[0][-1]
        if id not in self.data:
            self.data[id] = self._generate_data(rows[0][0])

        return self.data[id].copy()

    def get_multi(self, **options):
        """
        Gets row for sqlite3 database with specific variables specified in
        options. Same as get() except will return an list of all items that
        meet the selection criteria.

        :param options: sqlite3 select columns and values
        :return: unstringified data from database
        """
        query = self._generate_query(options)
        rows = self.connection.execute(query).fetchall()

        if len(rows) == 0:
            raise exceptions.SqliteLibraryException(
                "query returned no rows: " + self._args_to_str(options) )

        datas = []
        for r in rows:
            id = r[-1]
            if id not in self.data:
                self.data[id] = self._generate_data(r[0])
            datas.append(self.data[id].copy())
        return datas

    def get_random(self):
        """
        Gets a random row of data from sqlite3 database. Mostly used for
        testing purposes.

        :return: random row of unstringified data from database file
        """

        id = str(random.randint(1, self.max_size-1))
        return self.get(id=id)

    def load_all(self, limit=99999):
        """
        Loads all data from database into memory. This has to be done if one
        wants to iterate over all the rows in database. Can set a limit with
        limit=N, will only load that many elements.

        :param limit: the max number of database rows that should be loaded
        :type limt: int

        :return: None
        """

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
    """
    Sqlite3 library for holding stringifed motif.Motif objects. It is
    significantly faster then reloading each motif from file.

    :examples:

    ..  code-block:: python

        # gets all twoway junctions
        >>> mlib = sqlite_library.MotifSqliteLibrary("ideal_helices")
        >>> mlib.get(name="HELIX.IDEAL.2")
        <Motif(
	            structure='<Structure(name: N/A, #chains: 2, #residues: 8, #atoms: 172)>',
	            ends='2')>
    """

    def __init__(self, libname, rts=None):
        super(self.__class__, self).__init__(libname)
        self.rts = rts
        if self.rts is None:
            self.rts = residue_type.ResidueTypeSet()

    def get(self, **options):
        """
        Gets row for sqlite3 database with specific variables specified in
        options. Will at most return 1 item even if multiple items meet
        selection criteria.

        :param options: sqlite3 select columns and values
        :return: unstringified data from database
        """

        query = self._generate_query(options)
        rows = self.connection.execute(query).fetchall()

        if len(rows) == 0:
            raise exceptions.SqliteLibraryException(
                "query returned no rows: " + self._args_to_str(options) )

        id = rows[0][-1]
        if id not in self.data:
            self.data[id] = self._generate_data(rows[0][0])

        return  motif.Motif.copy(self.data[id], new_uuid=1)

    def get_multi(self, **options):
        """
        Gets row for sqlite3 database with specific variables specified in
        options. Same as get() except will return an list of all items that
        meet the selection criteria.

        :param options: sqlite3 select columns and values
        :return: unstringified data from database
        """
        query = self._generate_query(options)
        rows = self.connection.execute(query).fetchall()

        if len(rows) == 0:
            raise exceptions.SqliteLibraryException(
                "query returned no rows: " + self._args_to_str(options) )

        datas = []
        for r in rows:
            id = r[-1]
            if id not in self.data:
                self.data[id] = self._generate_data(r[0])
            datas.append(motif.Motif.copy(self.data[id], new_uuid=1))
        return datas

    def _generate_data(self, s):
        """
        Converts string to motif.Motif object

        :param s: String verision of motif.Motif object
        :type s: string

        :return: returns unstringified version of motif.Motif
        :rtype: motif.Motif
        """

        return motif.Motif.from_str(s, self.rts)

    @staticmethod
    def get_libnames():
        """
        get dictionary of valid registered library names / paths. Only library
        names registered in this corresponding function in children classes
        can be loaded.

        :return: Dictionary of names / pathes for each valid sqlite3 database
        :rtype: dict
        """

        libnames = {
            "ideal_helices"   : "/motif_libraries/ideal_helices.db",
            "ideal_helices_reversed" :  "/motif_libraries/ideal_helices_reversed.db",
            "twoway"          : "/motif_libraries/twoway.db",
            "tcontact"        : "/motif_libraries/tcontact.db",
            "hairpin"         : "/motif_libraries/hairpin.db",
            "nway"            : "/motif_libraries/nway.db",
            "unique_twoway"   : "/motif_libraries/unique_twoway.db",
            "bp_steps"        : "/motif_libraries/bp_steps.db",
            "helix"           : "/motif_libraries/helix.db",

        }

        return libnames


class MotifEnsembleSqliteLibrary(SqliteLibrary):
    """
    Sqlite3 library for holding stringifed motif_ensemble.MotifEnsemble objects.
    It is significantly faster then reloading each motif from file.

    :examples:
    """

    def __init__(self, libname, rts=None):
        super(self.__class__, self).__init__(libname)
        self.rts = rts
        if self.rts is None:
            self.rts = residue_type.ResidueTypeSet()

    def get(self, **options):
        """
        Gets row for sqlite3 database with specific variables specified in
        options. Will at most return 1 item even if multiple items meet
        selection criteria.

        :param options: sqlite3 select columns and values
        :return: unstringified data from database
        """

        query = self._generate_query(options)
        rows = self.connection.execute(query).fetchall()

        if len(rows) == 0:
            raise exceptions.SqliteLibraryException(
                "query returned no rows: " + self._args_to_str(options))

        id = rows[0][-1]
        if id not in self.data:
            self.data[id] = self._generate_data(rows[0][0])

        return motif_ensemble.MotifEnsemble.copy(self.data[id])

    def get_multi(self, **options):
        """
        Gets row for sqlite3 database with specific variables specified in
        options. Same as get() except will return an list of all items that
        meet the selection criteria.

        :param options: sqlite3 select columns and values
        :return: unstringified data from database
        """
        query = self._generate_query(options)
        rows = self.connection.execute(query).fetchall()

        if len(rows) == 0:
            raise exceptions.SqliteLibraryException(
                "query returned no rows: " + self._args_to_str(options))

        datas = []
        for r in rows:
            id = r[-1]
            if id not in self.data:
                self.data[id] = self._generate_data(r[0])
            datas.append(motif_ensemble.MotifEnsemble.copy(self.data[id]))
        return datas

    def _generate_data(self, s):
        return motif_ensemble.MotifEnsemble.from_str(s, self.rts)

    @staticmethod
    def get_libnames():
        libnames = {
            "bp_steps" :  "/motif_ensemble_libraries/bp_steps.db"
        }

        return libnames


class MotifStateSqliteLibrary(SqliteLibrary):

    def __init__(self, libname):
        super(self.__class__, self).__init__(libname)
        self._setup(libname)

    @staticmethod
    def get_libnames():
        libnames = {
            "ideal_helices" :  "/motif_state_libraries/ideal_helices.db",
            "twoway"        :  "/motif_state_libraries/twoway.db",
            "unique_twoway" :  "/motif_state_libraries/unique_twoway.db",
            "bp_steps"      :  "/motif_state_libraries/bp_steps.db",
            "twoway"        :  "/motif_state_libraries/twoway.db",
            "nway"          :  "/motif_state_libraries/nway.db",
            "hairpin"       :  "/motif_state_libraries/hairpin.db",
            "new_bp_steps"  :  "/motif_state_libraries/new_bp_steps.db",

        }

        return libnames

    def _generate_data(self, s):
        return motif.str_to_motif_state(s)

    def get(self, **options):
        m = super(self.__class__, self).get(**options)
        m.new_uuids()
        return m

    def to_motif_state_ensemble(self):
        self.load_all()
        mes = motif_ensemble.MotifStateEnsemble()
        motif_states = []
        energies = []
        for ms in self.all():
            motif_states.append(ms)
        mes.setup(self.name, motif_states, [1 for x in motif_states] )
        return mes


class MotifStateEnsembleSqliteLibrary(SqliteLibrary):
    def __init__(self, libname):
        super(self.__class__, self).__init__(libname)
        self._setup(libname)

    def _generate_data(self, s):
        return motif_ensemble.str_to_motif_state_ensemble(s)

    @staticmethod
    def get_libnames():
        libnames = {
            "bp_steps"      :  "/motif_state_ensemble_libraries/bp_steps.db",
            "twoway"        :  "/motif_state_ensemble_libraries/twoway.db",
        }

        return libnames


class MotifClusterSqliteLibrary(SqliteLibrary):

    def __init__(self, libname):
        super(self.__class__, self).__init__(libname)

    @staticmethod
    def get_libnames():
        libnames = {
            "twoways"      :  "/motif_state_libraries/twoways_clusters.db"
        }
        return libnames

    def _generate_data(self, s):
        return motif.str_to_motif_array(s)


def build_sqlite_library(path, data, keys, primary_key):
    """
    Builds a sqlite library to store various data. See setup scripts for
    more info on use. It is unlikely to be necessary to need to build new
    libraries. The length of each row must be the same length as keys.

    :param path: the path of sqlite3 library you want to create
    :param data: list of data to be stored in sqlite library
    :param keys: name of columns for each row in sqlite library
    :param primary_key: which of the keys specified will be the primary key,
        for sqlite3 index purposes, this must be a unique key, different for
        each row.

    :type path: String
    :type data: list
    :type keys: list
    :type primary_key: String

    :return: None

    :examples:

    ..  code-block:: python

        # build test library
        # each list is a row in sqlite db
        >>> data = [['the_word', 0], ['the', 1], ['hello', 2]]
        # the name of each column
        >>> keys = ['word', 'id']
        >>> sqlite_library.build_sqlite_library("test.db", data, keys, 'id')

        # read library
        >>> import sqlite3
        >>> conn = sqlite3.connect("test.db")
        >>> conn.execute("SELECT * from data_table WHERE word=\'the_word\'").fetchone()
        (u'the_word', u'0')
    """

    if len(data) == 0:
        raise exceptions.SqliteLibraryException("data must be longer then 0")

    if len(data[0]) != len(keys):
        raise exceptions.SqliteLibraryException(
            "length of each row must be the same length as keys")

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














