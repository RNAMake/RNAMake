import uuid


class BasepairType(object):
    WC = 0
    GU = 1
    NC = 2

class Basepair(object):
    __slots__= [
        "_uuid"
    ]

    def __init__(self, bp_uuid=None):
        self._uuid = bp_uuid
        if self._uuid is None:
            self._uuid = uuid.uuid1()

    def __eq__(self, other):
        return self._uuid == other._uuid

    def __ne__(self, other):
        return self._uuid != other._uuid


def calc_bp_name(res):
    res1, res2 = res

    res1_name = res1.chain_id+str(res1.num)+str(res1.i_code)
    res2_name = res2.chain_id+str(res2.num)+str(res2.i_code)

    if res1.chain_id < res2.chain_id:
        return res1_name+"-"+res2_name
    if res1.chain_id > res2.chain_id:
        return res2_name+"-"+res1_name

    if res1.num < res2.num:
        return res1_name+"-"+res2_name
    else:
        return res2_name+"-"+res1_name


