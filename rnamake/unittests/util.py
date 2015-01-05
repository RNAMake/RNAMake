import logging
from StringIO import StringIO

UnittestState = 0

class UnittestType(object):
    BASIC = 0
    ALL = 1

def get_log_output(func, args):
    logger = logging.getLogger()
    out = StringIO()
    stream_handler = logging.StreamHandler(out)
    logger.addHandler(stream_handler)
    func(args)
    output = out.getvalue().strip()
    return output


def supress_log_output(func, args):
    logging.disable(60)
    result = func(*args)
    logging.disable(0)
    return result
