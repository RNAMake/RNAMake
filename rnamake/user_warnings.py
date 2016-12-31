import warnings

__all__ = ['warn']

def warning_on_one_line(message, category, filename, lineno,
                        file=None, line=None):

    # original warnings do not have .class_name which are only user warnings
    # make sure to catch these and just use __name__ which is less useful
    try:
        s = '%s:%s: %s' %(category.class_name, category.__name__, message)
    except:
        s = '%s:%s: %s' %(category.__name__, category.__name__, message)

    return s


def warn(message, wtype):
    return warnings.warn(message, wtype)


warnings.formatwarning = warning_on_one_line


class RNAMakeWarning(UserWarning):
    class_name = __name__


class PDBFormatWarning(RNAMakeWarning):
    pass


class RNAStructureWarning(RNAMakeWarning):
    pass


class MotifWarning(RNAMakeWarning):
    pass


class MotifMergerWarning(RNAMakeWarning):
    pass