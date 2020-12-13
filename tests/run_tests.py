import pandas as pd
import subprocess
import os
import logging
import colorlog
import glob
import shutil
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="")
    optional_named = parser.add_argument_group('optional args')
    optional_named.add_argument('-t', '--test', help='perform specific test',
                        required=False, type=str,)


    args = parser.parse_args()
    return args


def init_logger(dunder_name) -> logging.Logger:
    log_format = (
        "[%(asctime)s " "%(name)s " "%(funcName)s] " "%(levelname)s " "%(message)s"
    )
    bold_seq = "\033[1m"
    colorlog_format = f"{bold_seq}" "%(log_color)s" f"{log_format}"
    cur_logger = logging.getLogger(dunder_name)
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

    cur_logger.addHandler(handler)
    cur_logger.setLevel(logging.INFO)
    return cur_logger


def get_cmd(fname):
    f = open(fname)
    lines = f.readlines()
    f.close()
    return lines[0]


def get_required_log_lines(fname):
    f = open(fname)
    required_log_lines = [line[9:].rstrip() for line in f.readlines()]
    f.close()
    return required_log_lines


def check_required_log_lines(name, required_log_lines, log_lines):
    seen = {}
    for line in log_lines:
        if line in required_log_lines:
            seen[line] = 1
    if len(seen) == len(required_log_lines):
        logger.info("test: {} PASSED!".format(name))
        return
    logger.warning("test: {} FAILED!".format(name))
    logger.warning("found {} out of {} required log statements".format(
            len(seen), len(required_log_lines)))
    logger.warning("missing the following log statements:")
    for line in required_log_lines:
        if line not in seen:
            logger.warning(line)


def get_expected_file_dict(fname):
    f = open(fname)
    expected_file_dict = {line.rstrip(): 1 for line in f.readlines()}
    f.close()
    return expected_file_dict


def get_produced_files():
    files = glob.glob("*")
    remove = [".", "..", "cmd", "inputs", "outputs", "org_outputs"]
    for r in remove:
        if r in files:
            files.remove(r)
    return files


def get_log_lines(cmd_str, row):
    try:
        output = subprocess.check_output(cmd_str, shell=True).decode("utf8")
    except subprocess.CalledProcessError:
        logger.error("test: {} did not run properly!".format(row["name"]))
        output = ""
    log_lines = [line[9:].rstrip() for line in output.split("\n")]
    return log_lines


logger = init_logger("RUN_TESTS")


def main():
    args = parse_args()
    cur_path = os.path.abspath(os.path.curdir)
    df = pd.read_csv("tests.csv")
    run_specific = False
    test_name = ""
    if args.test:
        run_specific = True
        test_name = args.test

    for i, row in df.iterrows():
        if run_specific:
            if test_name != row["name"]:
                continue
        os.chdir(row["name"])
        cmd_str = get_cmd("cmd/COMMAND")
        log_lines = get_log_lines(cmd_str, row)
        if len(log_lines) == 0:
            continue
        required_log_lines = get_required_log_lines("cmd/EXPECTED")
        check_required_log_lines(row["name"], required_log_lines, log_lines)
        files = get_produced_files()
        if not os.path.isdir("outputs"):
            os.mkdir("outputs")
        for f in files:
            shutil.move(f, "outputs/" + f)

        os.chdir(cur_path)


if __name__ == "__main__":
    main()
