

def update_motif_str(s):
    spl = s.split("&")
    path = spl[0]
    name = spl[1]
    score = spl[2]
    block_end_add = spl[3]
    mtype = spl[4]
    structure = spl[5]
    bp_strs = spl[6]
    end_indexes = spl[7]

    if len(spl) > 8:
        end_ids = spl[8]
        ss_str = spl[9]
        protein_beads = spl[10]

    new_bp_str = ""
    bp_spls = bp_strs.split("@")
    for bp_str in bp_spls:
        print bp_str

    new_s = ""