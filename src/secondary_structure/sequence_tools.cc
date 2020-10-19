//
// Created by Joseph Yesselman on 2019-04-11.
//

#include <secondary_structure/sequence_tools.h>
#include <secondary_structure/residue.h>

namespace secondary_structure {

void
get_res_types_from_sequence(
        String const & sequence,
        ResTypes & residue_types) {
    residue_types.resize(sequence.length());
    int i = 0;
    for(auto const & res_name : sequence) {
        residue_types[i] = convert_res_name_to_type(res_name);
        i++;
    }
}

int
find_gc_helix_stretches(
        PoseOP p,
        int length) {
    int count = 0;
    int stretches = 0;
    for(auto const & h : p->helices()) {
        count = find_longest_gc_helix_stretch_in_motif(h);
        if(count >= length) {
            stretches += 1;
            break;
        }
    }
    return stretches;
}

int
find_res_types_in_pose(
        PoseOP p,
        ResTypes const & residue_types) {
    int i = 0, j = 0, count = 0;
    for(auto const & c : p->chains()) {
        i = 0; j = 0;
        for(auto const & r : c->residues()) {
            if(r->res_type() == residue_types[j]) {
                j++;
                if(j == residue_types.size()) {
                    count += 1;
                }
            }
            else {
                j = 0;
            }
        }
    }
    return count;
}

int
find_res_types_in_motif(
        MotifOP m,
        ResTypes const & residue_types) {
    int i = 0, j = 0, count = 0;
    for(auto const & c : m->chains()) {
        i = 0; j = 0;
        for(auto const & r : c->residues()) {
            if(r->res_type() == residue_types[j]) {
                j++;
                if(j == residue_types.size()) {
                    count += 1;
                }
            }
            else {
                j = 0;
            }
        }
    }
    return count;
}


int
find_longest_gc_helix_stretch(
        PoseOP p) {
    int count = 0;
    int longest = 0;
    for(auto const & h : p->helices()) {
        count = find_longest_gc_helix_stretch_in_motif(h);
        if(count >= longest) {
            longest = count;
        }
    }
    return longest;

}

int
find_longest_gc_helix_stretch_in_motif(
        MotifOP m) {
    int count = 0;
    int longest = 0;
    for(auto const & r: m->chains()[0]->residues()) {
        if(r->res_type() == ResType::G || r->res_type() == ResType::C) {
            count += 1;
        }
        else {
            count = 0;
        }
        if(count >= longest) {
            longest = count;
        }
    }
    return longest;
}


ResType
get_complement_res_type(
        ResType r_type) {
    if(r_type == ResType::A) { return ResType::U; }
    else if(r_type == ResType::C) { return ResType::G; }
    else if(r_type == ResType::G) { return ResType::C; }
    else if(r_type == ResType::U) { return ResType::A; }
    else { throw secondary_structure::Exception("cannot get complement res type"); }

}

}

































