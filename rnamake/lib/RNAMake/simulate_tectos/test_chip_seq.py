import rnamake
import rnamake.prediction.motif_ensemble_tree as motif_ensemble_tree
import rnamake.prediction.motif_ensemble as motif_ensemble
import rnamake.motif_tree_state as motif_tree_state
import rnamake.secondary_structure_tree as secondary_structure_tree
import rnamake.motif as motif
import rnamake.util as util
import rnamake.motif_tree as motif_tree
import random
import copy
import math
import sys
import argparse

path = rnamake.settings.RESOURCES_PATH + "prediction/pdb_ensembles/"

def parse_args():
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-s', help='steps', required=False)
    parser.add_argument('-cseq', required=True)
    parser.add_argument('-css',required=True)
    parser.add_argument('-fseq', required=False)
    parser.add_argument('-fss', required=False)
    args = parser.parse_args()
    return args

def get_wc_flow():
    seq = 'CTAGGAATCTGGAAGTACCGAGGAAACTCGGTACTTCCTGTGTCCTAG'
    ss  = '((((((....((((((((((((....))))))))))))....))))))'

    return seq, ss


def mts_to_me(mts):
    me = motif_ensemble.MotifEnsemble()
    ms = motif_ensemble.MotifState(mts, 1.0)
    me.motif_states.append(ms)
    return me


def get_step_motifs(steps):
    motif_names = []
    for i in range(1,len(steps)):
        full_step = steps[i-1] + "=" + steps[i]
        motif_names.append(full_step)
    return motif_names


def get_ss_name(ss):
    name = ""
    for e in ss:
        if e == '(' or e == ')':
            name += 'P'
        elif e == '.':
            name += 'U'
        elif e == '&':
            name += '-'
    return name


def extract_steps_from_ss_tree(ss_tree):
    node = None
    for n in ss_tree.nodes:
        if n.ss_type == "Bulge":
            node = n
            break

    steps = []
    while node != None:
        if len(node.children) == 0:
            break
        c = node.children[0]
        if c.ss_type != "Basepair":
            break
        step = list(c.bp_type)
        if step[0] == "T":
            step[0] = "U"
        if step[1] == "T":
            step[1] = "U"
        steps.append("".join(step))
        node = c
    return steps

def replace_Ts_with_Us(seq):
    new_seq = ""
    for e in seq:
        if e == "T":
            new_seq += "U"
        else:
            new_seq += e
    return new_seq


def get_all_steps(ss_tree):
    node = None
    for n in ss_tree.nodes:
        if n.ss_type == "Bulge":
            node = n
            break

    required_nodes = []
    node = node.children[0]
    while len(node.children) > 0:
        required_nodes.append(node)
        node = node.children[0]

    i = 0
    steps = []
    while i < len(required_nodes)-1:
        if required_nodes[i+1].ss_type != "Basepair":
            bp1 = required_nodes[i].bp_type
            bp2 = required_nodes[i+2].bp_type
            motif_seq1 = ''.join(required_nodes[i+1].x_seq)
            motif_seq2 = ''.join(required_nodes[i+1].y_seq)

            seq = bp1[0] + motif_seq1 + bp2[0] + "-" + bp2[1] + motif_seq2[::-1] + bp1[1]
            seq = replace_Ts_with_Us(seq)
            ss = "("
            for j in range(len(motif_seq1)):
                ss += "."
            ss += '(&)'
            for j in range(len(motif_seq2)):
                ss += "."
            ss += ')'

            steps.append(path + seq + '_' +  get_ss_name(ss))
            i += 2

        else:
            step = required_nodes[i].bp_type + "=" + required_nodes[i+1].bp_type
            step = replace_Ts_with_Us (step)
            steps.append(step)

            i += 1

    return steps


def get_met(flow_ss_tree, chip_ss_tree):
    f = open('tetraloop.str')
    lines = f.readlines()
    f.close()

    ggaa_motif = motif.str_to_motif(lines[0])
    gaaa_motif = motif.str_to_motif(lines[1])

    ggaa_state = motif_tree_state.motif_to_state(ggaa_motif, end_index=1,
                                                end_flip=1)
    gaaa_state = motif_tree_state.motif_to_state(gaaa_motif, end_index=0,
                                                end_flip=1)

    flow_steps = extract_steps_from_ss_tree(flow_ss_tree)

    motif_names = get_step_motifs(flow_steps[1:])

    met = motif_ensemble_tree.MotifEnsembleTree()
    #flow
    met.add_ensemble(motif_ensemble.MotifEnsemble("AU=GC", 0, 0))
    met.add_ensemble(mts_to_me(ggaa_state))
    met.add_ensemble(motif_ensemble.MotifEnsemble(motif_names[0], 0, 1), parent_end_index=2)
    for i in range(1, len(motif_names)):
        met.add_ensemble(motif_ensemble.MotifEnsemble(motif_names[i], 0, 1))
    #chip
    met.add_ensemble(mts_to_me(gaaa_state))

    motif_names = get_all_steps(chip_ss_tree)[1:]
    for name in motif_names:
        print name
    print

    met.add_ensemble(motif_ensemble.MotifEnsemble(motif_names[0], 0, 1),
                      parent_end_index=1)
    for i in range(1, len(motif_names)):
        met.add_ensemble(motif_ensemble.MotifEnsemble(motif_names[i], 0, 1))
    return met


#mtst.last_node.states[1].d
def calculate_frame_score(target, target_flip, current):
    dist = util.distance(current.d, target.d)
    frame_score = dist
    r_diff = util.matrix_distance(current.r, target.r)
    r_diff_flip = util.matrix_distance(current.r, target_flip.r)
    if r_diff > r_diff_flip:
        r_diff = r_diff_flip
    frame_score += r_diff

    return frame_score


args = parse_args()
chip_ss_tree = secondary_structure_tree.SecondaryStructureTree(args.css, args.cseq)
if not args.fss:
    seq, ss = get_wc_flow()
    flow_ss_tree = secondary_structure_tree.SecondaryStructureTree(ss, seq)
else:
    flow_ss_tree = secondary_structure_tree.SecondaryStructureTree(args.fss,
                                                                   args.fseq)


steps = 100000
if args.s is not None:
    steps = int(args.s)

met = get_met(flow_ss_tree, chip_ss_tree)
mtst = met.get_mtst()
mtst.nodes_to_pdbs()
mtst.to_pdb('complex.pdb')
