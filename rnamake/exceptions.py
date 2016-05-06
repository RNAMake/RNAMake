
class GraphException(Exception):
    pass


class GraphIndexException(GraphException):
    def __init__(self, n_index, message):
        self.n_index = n_index
        super(self.__class__, self).__init__(message)


class GraphInvalidEndException(GraphException):
    def __init__(self, n, n_pos, message):
        self.n = n
        self.n_pos = n_pos
        super(self.__class__, self).__init__(message)


class TreeException(Exception):
    pass


class TreeIndexException(TreeException):
    def __init__(self, n_index, message):
        self.n_index = n_index
        super(self.__class__, self).__init__(message)

class TreeEndIndexException(TreeException):
    def __init__(self, n, n_index, message):
        self.n = n
        self.n_index = n_index
        super(self.__class__, self).__init__(message)

class SqliteLibraryException(Exception):
    pass

class MotifTypeException(Exception):
    pass

class X3dnaException(Exception):
    pass

class SecondaryStructureException(Exception):
    pass

class ResidueException(Exception):
    pass

class ChainException(Exception):
    pass

class StructureException(Exception):
    pass

class BasepairException(Exception):
    pass

class PDBParserException(Exception):
    pass

class MotifClustersException(Exception):
    pass