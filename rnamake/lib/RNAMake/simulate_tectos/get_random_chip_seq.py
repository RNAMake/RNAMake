import rnamake
import os
import random

def format_ss(ss_str):
    ss = ""
    halfway = 0
    for e in ss_str:
        if e == 'P' and not halfway:
            ss += '('
        if e == 'P' and halfway:
            ss += ')'
        if e == 'U':
            ss += '.'
        if e == '-':
            ss += '-'
            halfway = 1

    return ss

class SeqandSs(object):
    def __init__(self, seqs, ss):
        self.seqs, self.ss = seqs, ss

prime5_seq = 'CTAGGATATGG'
prime3_seq = 'CCTAAGTCCTAG'
prime5_ss  = '(((((((..(('
prime3_ss  = '))...)))))))'
tetra_seq  = 'GGGAAC'
tetra_ss   = '(....)'

ss = ['(', ')']
bps = [SeqandSs(['A','U'], ss),
       SeqandSs(['U','A'], ss),
       SeqandSs(['G','C'], ss),
       SeqandSs(['C','G'], ss)]


path = rnamake.settings.RESOURCES_PATH + "prediction/pdb_ensembles/"

twoways = []

for d in os.listdir(path):
    spl = d.split('.')
    if spl[-1] != "pop":
        continue

    name_spl = spl[0].split('_')
    seq_spl = name_spl[0].split('-')
    ss = format_ss(name_spl[1])
    ss_spl = ss.split('-')
    seqs = [seq_spl[0], seq_spl[1][::-1]]
    twoways.append(SeqandSs(seqs, ss_spl))

elements = []
for i in range(4):
    elements.append(random.choice(bps))
elements.append(random.choice(twoways))
for i in range(4):
    elements.append(random.choice(bps))

s1, s2, ss1, ss2 = "", "", "", ""
for e in elements:
    s1 += e.seqs[0]
    s2 += e.seqs[1]
    ss1 += e.ss[0]
    ss2 += e.ss[1]

s  = prime5_seq + s1 + 'GGGAAC' + s2[::-1] + prime3_seq
ss = prime5_ss  + ss1 +'(....)' + ss2  + prime3_ss

print s
print ss


