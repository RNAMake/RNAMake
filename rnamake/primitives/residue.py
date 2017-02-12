import abc
import uuid

import base

class Residue(base.BaseStructureObject):
    __metaclass__ = abc.ABCMeta

    __slots__ = [
        "_name",
        "_num",
        "_chain_id",
        "_i_code",
        "_uuid"]

    def __init__(self, name, num, chain_id, i_code=None, r_uuid=None):
        self._name = name
        self._num = num
        self._chain_id = chain_id
        self._i_code = i_code
        self._uuid = r_uuid

        if self._i_code is None:
            self._i_code = " "

        if self._uuid is None:
            self._uuid = uuid.uuid1()

    def __eq__(self, other):
        return self._uuid == other._uuid

    def __ne__(self, other):
        return self._uuid != other._uuid

    @property
    def name(self):
        return self._name

    @property
    def num(self):
        return self._num

    @property
    def chain_id(self):
        return self._chain_id

    @property
    def i_code(self):
        return self._i_code

    @property
    def uuid(self):
        return self._uuid