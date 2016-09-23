import sys
import random
import vienna_parameters
import os
import subprocess
import eterna_utils
import random
import re
from subprocess import Popen, PIPE, STDOUT
import thread, time, sys
from threading import Timer

import rnamake.settings

DEFAULT_TEMPERATURE = 37.0
BASES = ['A','U','G','C']
#fold a sequence
#@param seq:sequence
#@return: [parenthesis notation, energy]

def fold(seq):
    """
    folds sequence using Vienna

    args:
    seq is the sequence string

    returns:
    secondary structure
    """
    # run ViennaRNA
    if '&' in seq:
        p = Popen([rnamake.settings.VIENNA_BIN + 'RNAcofold', '-T','37.0'],
                  stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    else:
        p = Popen([rnamake.settings.VIENNA_BIN + 'RNAfold', '-T','37.0'],
                  stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    pair= p.communicate(input=''.join(seq))[0]
    p.wait()

    # split result by whitespace
    toks = re.split('\s+| \(?\s?',pair)
    ret= []
    ret.append(toks[1])
    ret.append(toks[2][1:-1])
    return ret

def fill_gc(elem , pair_map , seq, rand ):
    if(elem.type_ != eterna_utils.RNAELEMENT_STACK):
        return
    indices = elem.indices_
    length = len(indices)
    for ii in range(0,length):
        idx = indices[ii]
        if(pair_map[idx]<idx):
            continue
        if(rand.randint(0,1)==0):
            seq[idx]="G"
            seq[pair_map[idx]]="C"
        else:
            seq[idx]="C"
            seq[pair_map[idx]]="G"

def tout():
    thread.interrupt_main()

def timeout(func, args=(), timeout_duration=10, default=[]):
    ret = []
    finished=False
    seq=default
    t0=timeout
    try:
        timer = Timer(timeout_duration, tout)
        timer.start()
        t0=time.clock()
        seq = func(*args)
        finished=True
        t0=time.clock()-t0
        timer.cancel()
    except:
        print "time out"
    if(finished):
        ret.append(seq)
        ret.append(t0)
    else:
        ret.append(default)
        ret.append("timeout")
    return ret

