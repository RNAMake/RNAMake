import logging

def get_logger(name):
    logger = logging.getLogger(name)
    logger.propagate = False
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("[%(name)30s ] %(message)s")
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)
    return logger
