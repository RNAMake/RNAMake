import pandas as pd
import subprocess
import os
import logging
import colorlog

def init_logger(dunder_name) -> logging.Logger:
    log_format = (
        "[%(asctime)s " "%(name)s " "%(funcName)s] " "%(levelname)s " "%(message)s"
    )
    bold_seq = "\033[1m"
    colorlog_format = f"{bold_seq}" "%(log_color)s" f"{log_format}"
    logger = logging.getLogger(dunder_name)
    # colorlog.basicConfig(format=colorlog_format, datefmt="%H:%M")
    handler = colorlog.StreamHandler()
    handler.setFormatter(
        colorlog.ColoredFormatter(
            colorlog_format,
            datefmt="%H:%M",
            reset=True,
            log_colors={
                "DEBUG": "cyan",
                "WARNING": "yellow",
                "ERROR": "red",
                "CRITICAL": "red,bg_white",
            },
        )
    )

    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    return logger

def get_required_log_lines(fname):
    f = open(fname)
    required_log_lines = [l[9:].rstrip() for l in f.readlines()]
    f.close()
    return required_log_lines

def check_required_log_lines(name, required_log_lines, log_lines):
    seen = {}
    for l in log_lines:
        if l in required_log_lines:
            seen[l] = 1
    if len(seen) == len(required_log_lines):
        logger.info("test: {} PASSED!".format(name))
        return
    logger.warning("test: {} FAILED!".format(name))
    logger.warning("missing the following log statements:")
    for l in required_log_lines:
        if l not in seen:
            logger.warning(l)

logger = init_logger("RUN_TESTS")

def main():
    cur_path = os.path.abspath(os.path.curdir)
    df = pd.read_csv("tests.csv")
    for i, row in df.iterrows():
        os.chdir(row["name"])
        f = open("COMMAND")
        lines = f.readlines()
        f.close()
        required_log_lines = get_required_log_lines("EXPECTED")
        try:
            output = subprocess.check_output(lines[0], shell=True).decode("utf8")
        except:
            logger.error("test: {} did not run properly!".format(row['name']))
            continue
        log_lines = [l[9:].rstrip() for l in output.split("\n")]
        check_required_log_lines(row["name"], required_log_lines, log_lines)
        exit()


if __name__ == "__main__":
    main()
