#!/usr/bin/env python

import argparse
import os
import logging
import sys
from rnamake import resource_manager as rm
from rnamake import util, motif, motif_factory, motif_graph

FORMAT = "[%(asctime)s] %(message)s"
logging.basicConfig(format=FORMAT, stream=sys.stdout)
logger = logging.getLogger()
logger.setLevel(logging.INFO)

def parse_args():
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-pdb_50s', help='pdb of 50s ribosome', required=True)
    parser.add_argument('-pdb_30s', help='pdb of 30s ribosome', required=True)
    parser.add_argument('-cutpoint_50s',
                        help='which basepair to cut off a hairpin to build from',
                        default="A2854-A2863")
    parser.add_argument('-cutpoint_30s',
                        help= "which basepair to cut off a hairpin to build from",
                        default="A1445-A1457")
    args = parser.parse_args()
    return args


def build_motif_file(pdb_path, name):
    logger.info("building new motif file for: " + pdb_path + " will save to reduce run time "
                "for repeat runs. Can be overrided with -rebuild")

    rm.manager.add_motif(path=pdb_path, name=name, include_protein=1)
    m = rm.manager.get_motif(name=name)
    f = open("motifs/"+name+".motif", "w")
    f.write(m.to_str())
    f.close()

def load_motif_file(motif_file_path):
    logger.info("loading saved motif file: " + motif_file_path)
    m = motif.file_to_motif(motif_file_path)
    return m

def remove_hairpin_from_motif(m, bp_name):
    bps = m.get_basepair(name=bp_name)

    if len(bps) == 0:
        raise RuntimeError(
            "cannot find basepair with name: %s " % bp_name)

    bp = bps[0]
    bp.bp_type = "cW-W"

    removed = None
    for c in m.chains():
        try:
            removed = c.subchain(start_res=bp.res1, end_res=bp.res2)
        except:
            continue

    if removed is None:
        raise RuntimeError(
            "cannot remove a chain with both residues in the supplied basepair")

    res = removed.residues
    new_res = []
    new_bps = [bp]

    for r in m.residues():
        if r in bp.residues():
            new_res.append(r)
        if r in res:
            continue
        new_res.append(r)

    for bp in m.basepairs:
        if bp.res1 in res or bp.res2 in res:
            continue
        new_bps.append(bp)

    new_m = motif_factory.factory.motif_from_res(new_res, new_bps)
    new_m.name = m.name
    found = 0
    if bps[0] not in new_m.ends:
        new_m.ends.append(bps[0])
        motif_factory.factory._setup_secondary_structure(new_m)
    new_m.block_end_add = -1
    new_m = standardize_basepair_at_cutpoint(new_m)
    return new_m

def standardize_basepair_at_cutpoint(m):
    start_m = rm.manager.get_motif(name="HELIX.IDEAL.2")
    mg = motif_graph.MotifGraph()
    mg.add_motif(m)
    pos = mg.add_motif(start_m)
    if pos != -1:
        return m

    mg.get_node(0).data.ends[0].flip()
    pos = mg.add_motif(start_m)
    if pos != -1:
        return mg.get_node(0).data
    else:
        logger.warn("cannot standardize section of the ribosome: " + m.name + " this is "
                    "likely to lead to problem")


if __name__ == "__main__":
    args = parse_args()

    if not os.path.isdir("motifs"):
        logger.info("creating motifs/ directory to store motifs files")
        os.mkdir("motifs")

    name_50s = util.filename(args.pdb_50s)[:-4]
    name_30s = util.filename(args.pdb_30s)[:-4]
    if not os.path.isfile("motifs/"+name_50s+".motif"):
        build_motif_file(args.pdb_50s, name_50s)
    if not os.path.isfile("motifs/"+name_30s+".motif"):
        build_motif_file(args.pdb_30s, name_30s)

    m_50s = load_motif_file("motifs/"+name_50s+".motif")
    m_30s = load_motif_file("motifs/"+name_30s+".motif")

    m_50s_cut = remove_hairpin_from_motif(m_50s, "A2854-A2863")
    m_30s_cut = remove_hairpin_from_motif(m_30s, "A1445-A1457")

    mg = motif_graph.MotifGraph()
    mg.set_options('sterics', 0)
    mg.add_motif(m_30s_cut)
    mg.add_motif(m_50s_cut, orphan=1)
    mg.get_node(0).data.to_pdb("30S.pdb")
    mg.get_node(1).data.to_pdb("50S.pdb")

    f = open("test.mg", "w")
    f.write(mg.to_str() + "\n")
    f.write("1 " + "A2854-A2863" + "\n")
    f.write("0 " + "A1445-A1457" + "\n")
    f.close()

    logger.info("outputing start pdbs 30S.pdb and 50S.pdb")
    logger.info("success! COMMAND file generated. run \"source COMMAND\" to build tether")

    f = open("COMMAND", "w")
    f.write("design_rna -mg test.mg -only_ideal -design_pdbs -designs 1\n")
    f.close()




