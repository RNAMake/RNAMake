import abc

class BaseStructureObject(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def copy(cls, c):
        pass

    @abc.abstractmethod
    def from_str(cls, s, rts):
        pass

    @abc.abstractmethod
    def to_str(self):
        pass


class Transformable(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def move(self, p):
        pass

    @abc.abstractmethod
    def transform(self, t):
        pass
