import warnings

__all__ = ['warn']

def warning_on_one_line(message, category, filename, lineno,
                        file=None, line=None):
    return '%s:%s: %s' %(category.class_name, category.__name__, message)


def warn(message, wtype):
    return warnings.warn(message, wtype)


warnings.formatwarning = warning_on_one_line


class RNAMakeWarning(UserWarning):
    class_name = __name__

class PDBFormatWarning (RNAMakeWarning):
    pass