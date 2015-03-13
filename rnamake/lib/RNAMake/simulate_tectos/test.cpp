//
//  test.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 3/12/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//


#include <iostream>
#include <random>
#include <time.h>
#include "secondary_structure_tree.h"
#include "option.h"
#include "types.h"
#include "motif_ensemble.h"
#include "motif_ensemble_tree.h"
#include "motif_tree_state.h"
#include "motif_tree_state_tree.h"
#include "settings.h"

String mismatch_path = "/Users/josephyesselman/Downloads/results_1bp/";


Options
get_options(
    int argc,
    char const ** argv) {
    
    Strings allowed = split_str_by_delimiter("cseq,css,fseq,fss,steps", ",");
    return parse_command_into_options(argc, argv, allowed);
    
    
}

StringStringMap
get_construct_seq_and_ss(
     Options const & options) {
    
    //defaults are wc flow and wc chip (id: 1)
    StringStringMap defaults;
    defaults["fseq"] = "CTAGGAATCTGGAAGTACCGAGGAAACTCGGTACTTCCTGTGTCCTAG";
    defaults["fss" ] = "((((((....((((((((((((....))))))))))))....))))))";
    defaults["cseq"] = "CTAGGATATGGAAGATCCTCGGGAACGAGGATCTTCCTAAGTCCTAG";
    defaults["css" ] = "(((((((..((((((((((((....))))))))))))...)))))))";
    Strings keys = split_str_by_delimiter("fseq,fss,cseq,css", ",");
    StringStringMap constructs;
    int found = 0;
    for (auto const & k : keys) {
        found = 0;
        for (auto const & o : options) {
            if(o.key.compare(k) == 0) {
                found = 1;
                constructs[k] = o.value;
                break;
            }
        }
        
        if(!found) { constructs[k] = defaults[k]; }
    }
    return constructs;
    
}

String
remove_Ts(String const & step) {
    String nstep;
    for(auto const & s : step) {
        if(s == 'T') { nstep += 'U'; }
        else         { nstep += s;   }
    }
    return nstep;
}

Strings
get_steps_from_ss_tree(
    SecondaryStructureTree const & ss_tree) {
    SecondaryStructureNodeOPs bulges = ss_tree.get_bulges();
    SecondaryStructureNodeOP current = bulges[0]->children()[0];
    Strings steps;
    SecondaryStructureNodeOPs required_nodes;
    int bulge_flank_size = 2;
    while (current->ss_type() != SSN_HAIRPIN) {
        required_nodes.push_back(current);
        current = current->children()[0];
    }
    int i = 1;
    String step;
    while(i < required_nodes.size()) {
        if(i+bulge_flank_size-1 < required_nodes.size() &&
           required_nodes[i+bulge_flank_size-1]->ss_type() == SSN_TWOWAY) {
            String first_bp = required_nodes[i-1]->bp_type();
            String second_bp = required_nodes[i]->bp_type();
            String third_bp = required_nodes[i+2]->bp_type();
            String forth_bp = required_nodes[i+3]->bp_type();
            Strings seqs = split_str_by_delimiter(required_nodes[i+1]->seq(), "+");
            String seq1, seq2;
            seq1.push_back(second_bp[0]);
            seq1 += seqs[0];
            seq1.push_back(third_bp[0]);
            seq2.push_back(second_bp[1]);
            seq2 += seqs[1];
            seq2.push_back(third_bp[1]);
            std::reverse(seq2.begin(), seq2.end());
            String name = seq1 + "-" + seq2;
            std::cout << name << std::endl;
            exit(0);
            steps.push_back(remove_Ts(first_bp+"="+second_bp));
            //std::cout << name << std::endl;
            steps.push_back(remove_Ts(name));
            steps.push_back(remove_Ts(third_bp+"="+forth_bp));
            i = i+4;
            continue;
        }
        step = required_nodes[i-1]->bp_type() + "=" + required_nodes[i]->bp_type();
        steps.push_back(remove_Ts(step));
        i++;
        
    }
    //first one is in tetraloop receptor structure
    steps.erase(steps.begin());
    
    return steps;
}

Strings
get_lines_from_file(String const fname) {
    String line;
    Strings lines;
    std::ifstream input;
    input.open(fname);
    while ( input.good() ) {
        getline(input, line);
        if( line.length() < 10 ) { break; }
        lines.push_back(line);
        
    }
    return lines;
}

MotifEnsemble
mts_to_me(MotifTreeStateOP const & mts) {
    MotifEnsemble me;
    MotifState ms (mts, 1.0);
    me.add_motif_state(ms);
    return me;
}

MotifEnsemble
get_ensemble(
             String const & name,
             int end,
             int flip) {
    if (name[2] == '=') {
        return MotifEnsemble(name, end, flip);
    }
    else {
        String path = mismatch_path + name;
        return MotifEnsemble(path, end, flip);
    }
}


MotifEnsembleTree
get_met(
        SecondaryStructureTree const & flow_ss_tree,
        SecondaryStructureTree const & chip_ss_tree) {
    
    Strings lines = get_lines_from_file(base_dir() + "/rnamake/lib/RNAMake/simulate_tectos/tetraloop.str");
    ResidueTypeSet rts;
    MotifOP ggaa_motif ( new Motif(lines[0], rts));
    MotifOP gaaa_motif ( new Motif(lines[1], rts));
    MotifTreeStateOP ggaa_state = motif_to_state(ggaa_motif, 1, 1);
    MotifTreeStateOP gaaa_state = motif_to_state(gaaa_motif, 0, 1);
    Strings flow_steps = get_steps_from_ss_tree(flow_ss_tree);
    Strings chip_steps = get_steps_from_ss_tree(chip_ss_tree);
    MotifEnsembleTree met;
    met.add_ensemble(MotifEnsemble("AU=GC", 0, 0));
    met.add_ensemble(mts_to_me(ggaa_state));
    met.add_ensemble(MotifEnsemble(flow_steps[0], 0, 1), NULL, 2);
    for (int i = 1; i < flow_steps.size(); i++) {
        met.add_ensemble(MotifEnsemble(flow_steps[i], 0, 1));
    }
    met.add_ensemble(mts_to_me(gaaa_state));
    met.add_ensemble(get_ensemble(chip_steps[0], 0, 1), NULL, 1);
    for (int i = 1; i < chip_steps.size(); i++) {
        met.add_ensemble(get_ensemble(chip_steps[i], 0, 1));
    }
    
    return met;
}


int main(int argc, const char * argv[]) {
    Options options = get_options(argc, argv);
    StringStringMap constructs_seq_and_ss = get_construct_seq_and_ss(options);
    constructs_seq_and_ss["cseq"]= "CTAGGATATGGGGCCGAUGGUGGGAACACCAUACCGGCCCCTAAGTCCTAG";
    constructs_seq_and_ss["css" ]= "(((((((..(((((((((((((....))))))..)))))))...)))))))";
    SecondaryStructureTree chip_ss_tree ( constructs_seq_and_ss["css"], constructs_seq_and_ss["cseq"]);
    SecondaryStructureTree flow_ss_tree ( constructs_seq_and_ss["fss"], constructs_seq_and_ss["fseq"]);
    MotifEnsembleTree met = get_met(flow_ss_tree, chip_ss_tree);
    int nsteps = 10000000;
    return 0;
}


