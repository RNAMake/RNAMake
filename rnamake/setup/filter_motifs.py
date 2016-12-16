import os
import glob
import pandas as pd
from rnamake import pose_factory, motif_type, motif_factory, util, motif

def build_structure_motifs():

    path = "/Users/josephyesselman/projects/REDESIGN/resources/non-redundant-rnas/"
    dirs = list(os.walk(path))[0][1]
    for d in dirs:
        try:
            m = motif_factory.factory.motif_from_file(path+d, include_protein=1)
        except:
            continue

        print m.name, len(m.protein_beads)
        f = open("structure_motifs/"+d+".motif", "w")
        f.write(m.to_str())
        f.close()

def build_all_motifs():
    path = "/Users/josephyesselman/projects/REDESIGN/resources/non-redundant-rnas/"

    types = [motif_type.HELIX, motif_type.TWOWAY, motif_type.NWAY, motif_type.HAIRPIN]

    motif_files = glob.glob("structure_motifs/*")
    for f in motif_files:
        name = util.filename(f)[:-6]
        print name
        p = pose_factory.factory.pose_from_file(path+name)
        for t in types:
            motifs = p.motifs(t)
            if len(motifs) == 0:
                continue
            t_name = motif_type.type_to_str(t)
            m_path = "sectioned_motifs/"+t_name+"/"+name+".motif"
            f = open(m_path, "w")
            for i, m in enumerate(motifs):
                m.name = t_name + "." + name + "." + str(i)
                f.write(m.to_str()+"\n")
            f.close()

def analyze_motifs():
    twoway_files = glob.glob("sectioned_motifs/TWOWAY/*")

    df = pd.DataFrame(columns = "m_name m_size extra_bps bps_names".split())
    pos = 0
    for m_file in twoway_files:
        f = open(m_file)
        lines = f.readlines()
        f.close()
        name = util.filename(m_file)[:-6]

        full_m = motif.file_to_motif("structure_motifs/"+name+".motif")

        for l in lines:
            m = motif.str_to_motif(l)
            res = m.residues()
            centers = [ util.center(r.atoms) for r in res]
            keys = { r.chain_id + " " + r.name + " " + str(r.num) : r for r in res }

            bps_names = ""
            found = 0
            p_clashes = 0
            for bp in full_m.basepairs:
                r1_key = bp.res1.chain_id + " " + bp.res1.name + " " + str(bp.res1.num)
                r2_key = bp.res2.chain_id + " " + bp.res2.name + " " + str(bp.res2.num)
                if r1_key in keys and r2_key not in keys:
                    bps_names += bp.name() + "|"
                    found += 1

                if r2_key in keys and r1_key not in keys:
                    bps_names += bp.name() + "|"
                    found += 1

            for b in full_m.protein_beads:
                for c in centers:
                    dist = util.distance(b.center, c)
                    if dist < 6:
                        p_clashes += 1

            if p_clashes > 1:
                print m.name

            m.to_pdb("twoways/"+m.name+".pdb")
            df.loc[pos] = [m.name, len(m.residues()), found, bps_names]
            pos += 1

    df.to_csv("motifs_extra_bps.csv", index=False)


def compile_summary():
    df = pd.read_csv("motifs_extra_bps.csv")

    #print len(df[df.apply(lambda x: x['extra_bps'] == 0, axis=1)])
    #print len(df[df.apply(lambda x: x['extra_bps'] == 1, axis=1)])
    #print len(df[df.apply(lambda x: x['extra_bps'] == 2, axis=1)])
    #print len(df[df.apply(lambda x: x['extra_bps'] > 2, axis=1)])

    print df[df.apply(lambda x: x['extra_bps'] == 0 and x['m_size'] > 15, axis=1)]

#need to check for more than 2 chains


#build_structure_motifs()
#build_all_motifs()
analyze_motifs()
#compile_summary()


"""if d != "1S72":
     continue
p = pose_factory.factory.pose_from_file(path+d)
two_ways_motifs = p.motifs(motif_type.TWOWAY)"""

