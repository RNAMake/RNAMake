//
//  main.cpp
//  motif_tree_path_refiner
//
//  Created by Joseph Yesselman on 3/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include <random>
#include <fstream>
#include "motif_library.h"
#include "motif_tree_state_selector.h"
#include "motif_tree_state_search_scorer.h"
#include "motif_tree_state_search_node.h"
#include "motif_tree_state_library.h"
#include "motif_tree_state_tree.h"
#include "motif_tree_state_search.h"
#include "motif_tree_state_path_refiner.h"
#include "secondary_structure_tree.h"
#include "motif_ensemble_tree.h"
#include "option.h"
#include "mtst_to_met_converter.h"
#include "sequence_designer.h"
#include "FileIO.h"
#include "vienna.h"
#include "motif_tree_state_alt_pather.h"

String mismatch_path = "/Users/josephyesselman/Downloads/results_1bp/";

MotifTreeStateTree
get_two_way_helix_mts_tree(int size=10) {;
    std::vector<MotifTreeStateLibrary> mts_libs(2);
    mts_libs[0] = MotifTreeStateLibrary ( HELIX );
    mts_libs[1] = MotifTreeStateLibrary ( TWOWAY );
    MotifTreeStateTree mtst;
    
    srand(unsigned(time(NULL)));
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<float> dist(0,1);
    
    int i = 0;
    int pos = 0;
    int count = 0;
    while (i < size) {
        if( i % 2 == 0) { pos = 0; }
        else            { pos = 1; }
        
        int mts_pos = dist(mt)*(mts_libs[pos].motif_tree_states().size()-1);
        MotifTreeStateOP mts = mts_libs[pos].motif_tree_states()[mts_pos];
        
        MotifTreeStateNodeOP node = mtst.add_state(mts, NULL);
        if (node != NULL) { i++; }
        count++;
        if(count > 10000) { break; }
        
    }
    
    return mtst;
    
}

MotifTreeStateTree
get_set_problem(Strings const & names) {
    
    std::vector<MotifTreeStateLibrary> mts_libs(2);
    mts_libs[0] = MotifTreeStateLibrary ( HELIX );
    mts_libs[1] = MotifTreeStateLibrary ( TWOWAY );
    MotifTreeStateTree mtst;
    
    int i = 0;
    int pos = 0;
    for(auto const & name : names) {
        if( i % 2 == 0) { pos = 0; }
        else            { pos = 1; }

        MotifTreeStateOP mts = mts_libs[pos].get_state(name);
        MotifTreeStateNodeOP node = mtst.add_state(mts, NULL);
        i++;
    }
    
    return mtst;
    
}

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
            seq1.push_back(first_bp[0]); seq1.push_back(second_bp[0]);
            seq1 += seqs[0];
            seq1.push_back(third_bp[0]); seq1.push_back(forth_bp[0]);
            seq2.push_back(first_bp[1]); seq2.push_back(second_bp[1]);
            seq2 += seqs[1];
            seq2.push_back(third_bp[1]); seq2.push_back(forth_bp[1]);
            std::reverse(seq2.begin(), seq2.end());
            String name = seq1 + "-" + seq2;
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
        if( line.length() < 2 ) { break; }
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
    return met;
}


int
test_path_refiner() {
    MotifTreeStateTree mtst = get_two_way_helix_mts_tree(10);
    mtst.to_pdb("test.pdb");
    MotifTreeStatePathRefiner path_finder;
    path_finder.set_numeric_option("max_size", mtst.last_node()->size());
    MotifTreeStateSearchSolutionOPs solutions = path_finder.find_path(mtst.nodes()[0]->states()[0], mtst.nodes().back()->active_states()[0]);
    
    return 0;
}

int
test_path_refiner_ttr() {
    
    Strings lines = get_lines_from_file(base_dir() + "/rnamake/lib/RNAMake/simulate_tectos/tetraloop.str");
    ResidueTypeSet rts;
    MotifOP gaaa_motif ( new Motif(lines[1], rts));
    Beads beads = gaaa_motif->get_beads(gaaa_motif->ends());
    Points centers;
    for (auto const & b : beads) {
        if(b.btype() == PHOS) { continue; }
        centers.push_back(b.center());
    }
    
    gaaa_motif->to_pdb("test.pdb");
    BasepairStateOP start (new BasepairState(gaaa_motif->ends()[1]->state()));
    BasepairStateOP end (new BasepairState(gaaa_motif->ends()[0]->state()));

    
    MotifTreeStateSearch search;
    MotifTreeStateLibraryOPs mts_libs(2);
    mts_libs[0] = MotifTreeStateLibraryOP ( new MotifTreeStateLibrary(HELIX, 5));
    mts_libs[1] = MotifTreeStateLibraryOP ( new MotifTreeStateLibrary(TWOWAY));
    MotifTreeStateSelectorOP selector ( new MotifTreeStateSelector(mts_libs, "round_robin"));
    search.base_beads(centers);
    search.set_numeric_option("max_size", 130);
    search.set_numeric_option("accept_ss_score", 0);
    search.set_numeric_option("max_n_solutions", 100000000);
    search.set_numeric_option("log", 1);
    search.set_numeric_option("save_solutions", 0);
    search.set_numeric_option("accept_score", 10);
    //search.set_numeric_option("max_node_level", 5);

    MotifTreeStateSearchScorerOP scorer ( new MTSS_Astar() );
    MotifTreeStateSearchSolutionOPs solutions = search.search(start, end, selector, scorer);
    
    std::stringstream ss;
    
    std::cout << solutions.size() << std::endl;
    
    int i = 0;
    for (auto const & s : solutions) {
        for (auto const & n : s->path() ) {
            std::cout << n->mts()->name() << " ";
        }
        std::cout << std::endl;
        ss << "solutions." << i << ".pdb";
        s->to_mtst().to_pdb(ss.str());
        ss.str("");
        i++;
    }
    
    return 1;
}

int
test_path_refiner_tecto(int argc, const char * argv[]) {
    Options options = get_options(argc, argv);
    StringStringMap constructs_seq_and_ss = get_construct_seq_and_ss(options);
    //constructs_seq_and_ss["cseq"]= "CTAGGATATGGAAGATACTCGGGAACGAGAATCTTCCTAAGTCCTAG";
    //constructs_seq_and_ss["css" ]= "(((((((..(((((((.((((....)))).)))))))...)))))))";
    SecondaryStructureTree chip_ss_tree ( constructs_seq_and_ss["css"], constructs_seq_and_ss["cseq"]);
    SecondaryStructureTree flow_ss_tree ( constructs_seq_and_ss["fss"], constructs_seq_and_ss["fseq"]);
    MotifEnsembleTree met = get_met(flow_ss_tree, chip_ss_tree);
    
    MotifTreeStateTree mtst = met.get_mtst();
    Points beads;
    int i = -1;
    for( auto const & n : mtst.nodes()) {
        i++;
        if( i == 1) { continue; }
        if(n == mtst.last_node()) { continue; }
        for( auto const & b : n->beads() ) {
            beads.push_back(b);
        }
    }
    mtst.to_pdb("start.pdb");
    BasepairStateOP target = mtst.nodes()[2]->states()[0];
    MotifTreeStateSearch search;
    //search.base_beads(beads);
    search.set_numeric_option("accept_score", 5.0);
    search.set_numeric_option("max_n_solutions", 1000);
    search.set_numeric_option("max_size", 60);
    MotifTreeStateSearchScorerOP scorer ( new MTSS_Astar() );
    MotifTreeStateSearchSolutionOPs solutions = search.search(mtst.last_node()->states()[1], target, NULL, scorer);
    std::cout << solutions.size() << std::endl;
    i = 0;
    std::stringstream ss;
    for(auto const & s : solutions) {
        std::cout << i << " | ";
        for(auto const & n : s->path()) { std::cout << n->mts()->name() << " "; }
        std::cout << " | " << s->score() << std::endl;
        ss << "solution." << i << ".pdb";
        s->to_mtst().to_pdb(ss.str());
        ss.str("");
        i++;
    }
    
    return 1;
}

int
test_converter() {
    //Strings mts_names = split_str_by_delimiter("HELIX.LE.7-0-0-0-0-1-0 TWOWAY.1S72.6-0-0-1-0-0-1 HELIX.LE.5-0-0-0-0-1-0 TWOWAY.4P95.0-0-0-1-0-0-1 HELIX.LE.20-0-0-0-0-1-1 TWOWAY.2VQE.35-0-0-1-0-0-0 HELIX.LE.12-0-0-0-0-1-1 TWOWAY.2HW8.0-0-0-1-0-0-0 HELIX.LE.15-0-0-0-0-1-1 TWOWAY.2VQE.26-0-0-0-0-1-1", " ");
    MotifTreeStateTree mtst = get_two_way_helix_mts_tree(10);
    for(auto const & n : mtst.nodes()) {
        std::cout << n->mts()->name() << std::endl;
    }
    //MotifTreeStateTree mtst = get_set_problem(mts_names);
    //mtst.to_pdb("test.pdb");
    /*MTSTtoMETConverter converter;
    MotifEnsembleTreeOP met = converter.convert(mtst);
    MotifTreeStateTree mtst2 = met->get_mtst();
    MotifTree mt = mtst2.to_motiftree();
    */
    /*for(auto const & n : mtst.nodes()) {
        std::cout << n->index() << " " << n->mts()->name() << std::endl;
    }
    
    std::cout << std::endl;
    for(auto const & n : mtst2.nodes()) {
        std::cout << n->index() << " " << n->mts()->name() << std::endl;
    }*/
    
    //mtst2.to_pdb("test2.pdb");
    
    return 1;
}

int
test_sequence_designer() {
    Strings lines = get_lines_from_file("mtss_ttr.log");
    Strings names;
    String ss, dseq, seq;
    float score;
    std::vector<MotifTreeStateLibrary> mts_libs(2);
    mts_libs[0] = MotifTreeStateLibrary(HELIX);
    mts_libs[1] = MotifTreeStateLibrary(TWOWAY);
    Strings tetraloops = get_lines_from_file(base_dir() + "/rnamake/lib/RNAMake/simulate_tectos/tetraloop.str");
    ResidueTypeSet rts;
    MotifOP gaaa_tetraloop ( new Motif(tetraloops[1], rts));
    int j = -1;
    
    MotifTreeStateOP gaaa_mts = motif_to_state(gaaa_tetraloop, 2);
    MotifTreeStateTree mtst;
    mtst.add_state(mts_libs[0].get_state("HELIX.LE.3-0-0-0-0-1-1"), NULL);
    mtst.add_state(gaaa_mts, NULL);
    mtst.sterics(0);
    std::vector<float> scores;

    SequenceDesigner designer;

    for( auto const & l : lines) {
        j++;
        int i = 0, pos = 0;
        
        i = -1;
        names = split_str_by_delimiter(l, " ");
        if(names.size() < 4) { continue; }
        for(auto const & name : names) {
            i++;
            if ( i == 0 ) { continue; }
            if ( i % 2 == 0)  { pos = 1; }
            else              { pos = 0; }
            
            if(i == 1) {
                mtst.add_state(mts_libs[pos].get_state(name), NULL, 1);
            }
            else {
                mtst.add_state(mts_libs[pos].get_state(name), NULL);
            }
            
        }
        
        
        //mtst.to_pdb("test.pdb");
        MotifTree mt = mtst.to_motiftree();
        
        while(mtst.last_node()->index() > 2) { mtst.remove_node(); }
        
        mt._add_connection(mt.nodes()[2], mt.nodes().back(), 1000);
        PoseOP p = mt.to_pose();
        if(p->chains().size() > 1) { continue; }
        ss   = p->secondary_structure();
        dseq = p->designable_sequence();
        seq = designer.design(dseq, ss);
        scores = designer.scores();
        std::cout << seq << std::endl;
        std::cout << j << " " << designer.score() << " ";
        for( auto const & s : scores ) { std::cout << s << " "; }
        std::cout << std::endl;
        exit(0);

        //ss = v.get_structure();
        

        //sequences.push_back(seq);
    }
    
    return 0;
}




int optimize(int j_pos) {
    
    Strings lines = get_lines_from_file("mtss_ttr.log");
    Strings names;
    String ss, dseq, seq;
    float score;
    std::vector<MotifTreeStateLibrary> mts_libs(2);
    mts_libs[0] = MotifTreeStateLibrary(HELIX);
    mts_libs[1] = MotifTreeStateLibrary(TWOWAY);
    Strings tetraloops = get_lines_from_file(base_dir() + "/rnamake/lib/RNAMake/simulate_tectos/tetraloop.str");
    ResidueTypeSet rts;
    MotifOP gaaa_tetraloop ( new Motif(tetraloops[1], rts));
    int j = -1;
    
    MotifTreeStateOP gaaa_mts = motif_to_state(gaaa_tetraloop, 2);
    MotifTreeStateTree mtst;
    mtst.add_state(mts_libs[0].get_state("HELIX.LE.3-0-0-0-0-1-1"), NULL);
    mtst.add_state(gaaa_mts, NULL);
    mtst.sterics(0);
    std::vector<float> scores;
    
    SequenceDesigner designer;
    MotifTreeStateAltPather alt_pather;
    std::vector<MotifTreeStateTree> all_trees;
    SeqandScores top_seqs;
    
    String prime5_seq = "GGAACAGCUCGAGUAGAGCUGAAA";
    String prime5_ss  = "....((((((.....))))))...";
    
    String prime3_seq = "AAACGCAUCGAGUAGAUGCGAACAAAGAAACAACAACAACAAC";
    String prime3_ss  = "...((((((.....)))))).......................";
    
    String common_seq_1 = "GGGAAAAGAGGAAAGGCGAAGGAAG";
    String common_seq_2 = "GGGAAUAGGGCAAGGGAAUAGGAGG";
    String common_ss    = ".........................";
    
    String submission_seq, submission_ss;
    String full_seq, full_ss;
    String designed_seq;

    int size_diff;
    int common_seq_size;
    int actual_start;
    
    MTSTtoMETConverter converter;
    MotifTreeStateNodeOP last;
    float dist;
    int count = 0;

    
    for( auto const & l : lines) {
        j++;
        int i = 0, pos = 0;
        
        i = -1;
        names = split_str_by_delimiter(l, " ");
        if(names.size() < 4) { continue; }
        for(auto const & name : names) {
            i++;
            if ( i == 0 ) { continue; }
            if ( i % 2 == 0)  { pos = 1; }
            else              { pos = 0; }
            
            if(i == 1) {
                mtst.add_state(mts_libs[pos].get_state(name), NULL, 1);
            }
            else {
                mtst.add_state(mts_libs[pos].get_state(name), NULL);
            }
            
        }
        
        
        MotifTree mt = mtst.to_motiftree();
        
        
        mt._add_connection(mt.nodes()[2], mt.nodes().back(), 1000);
        PoseOP p = mt.to_pose();
        if(p->chains().size() > 1) {
            while(mtst.last_node()->index() > 2) { mtst.remove_node(); }
            continue;
        
        }
        
        //all_trees = alt_pather.get_alt_paths(mtst);
        
        ss   = p->secondary_structure();
        dseq = p->designable_sequence();
        
        if(dseq.size() > 120) {
            while(mtst.last_node()->index() > 2) { mtst.remove_node(); }
            continue;
        }
        
        size_diff = 150 - (int)dseq.length();
        common_seq_size = size_diff / 2;
        submission_seq = common_seq_1.substr(0, common_seq_size) + dseq;
        submission_ss = common_ss.substr(0, common_seq_size) + ss;
        size_diff = 150 - (int)submission_seq.length();
        submission_seq = submission_seq + common_seq_2.substr(common_seq_2.size() - size_diff);
        submission_ss = submission_ss + common_ss.substr(0, size_diff);

        actual_start = common_seq_size + (int)prime5_seq.length();
        
        full_seq = prime5_seq + submission_seq + prime3_seq;
        full_ss  = prime5_ss  + submission_ss  + prime3_ss;
        
        designer.design(full_seq, full_ss);
        top_seqs = designer.top_seqs();
     
        
        //std::cout << ss << std::endl;
        //std::cout << top_seqs[0].seq << std::endl;
        //std::cout << top_seqs[0].score << std::endl;
   
        for( int k = 0; k < 10; k++) {
        
            designed_seq = top_seqs[k].seq.substr(actual_start, dseq.length());
            MotifEnsembleTreeOP met = converter.convert(mtst,mt,designed_seq);
            MotifTreeStateTree mtst2 = met->get_mtst();
            dist = mtst2.nodes()[4]->active_states()[0]->d().distance(mtst2.last_node()->active_states()[0]->d());
            std::cout << j << " " << top_seqs[k].seq.substr(24, 150) << " " << designed_seq << " " << top_seqs[k].score << " " <<  dist << std::endl;
        
        //mtst.to_pdb("org.pdb");
        //mtst2.to_pdb("converted.pdb");
        }
        
        while(mtst.last_node()->index() > 2) { mtst.remove_node(); }
        
        exit(0);

        
    }
    

    
    return 1;
}

MotifTreeStateSelectorOP
custom_selector(
    String lib1) {
    
    MotifTreeStateLibraryOPs mts_libs(2);
    mts_libs[0] = MotifTreeStateLibraryOP ( new MotifTreeStateLibrary(TWOWAY));
    mts_libs[1] = MotifTreeStateLibraryOP ( new MotifTreeStateLibrary(lib1));
    MotifTreeStateLibraryOP hlib ( new MotifTreeStateLibrary (HELIX, 5));
    
    SelectorNodeOPs nodes;
    SelectorNodeOP node (new SelectorNode(hlib, 0));
    nodes.push_back(node);
    int i = 1;
    for( auto const & mts_lib : mts_libs) {
        SelectorNodeOP n (new SelectorNode(mts_lib, i));
        nodes[0]->add_connection(n);
        n->add_connection(nodes[0]);
        nodes.push_back(n);
        i++;
    }
    
    nodes[2]->max_uses(1);
    nodes[2]->required_uses(1);

    MotifTreeStateSelectorOP selector ( new MotifTreeStateSelector(nodes));
    return selector;

    
}


int build_ribozyme() {
    
    Strings motifs = get_lines_from_file(base_dir() + "/rnamake/lib/RNAMake/motif_tree_path_refiner/hammerhead_motifs.str");
    ResidueTypeSet rts;
    MotifOP hammerhead_tc ( new Motif(motifs[0], rts));
    MotifOP hammerhead_core ( new Motif(motifs[1], rts));
    MotifOP fmn (new Motif(motifs[2], rts));

    MotifTreeStateLibraryOPs mts_libs(2);
    mts_libs[0] = MotifTreeStateLibraryOP ( new MotifTreeStateLibrary(HELIX, 5));
    mts_libs[1] = MotifTreeStateLibraryOP ( new MotifTreeStateLibrary(TWOWAY));
    
    MotifTree mt;
    MotifLibrary mlib ( HELIX );
    MotifTreeNodeOP node = mt.add_motif(hammerhead_core, NULL, 2);
    mt.add_motif(mlib.get_motif("HELIX.LE.20"));
    mt.add_motif(hammerhead_tc, NULL, 1);
    mt.add_motif(fmn, node);

    PoseOP start_conf = mt.to_pose();

    start_conf->to_pdb("test.pdb");
    Beads beads = start_conf->get_beads(start_conf->ends());
    Points centers;
    for (auto const & b : beads) {
        if(b.btype() == PHOS) { continue; }
        centers.push_back(b.center());
    }

    int i = 0;
    for (auto const & end : start_conf->ends()) {
        std::cout << i << " " << end->name() << std::endl;
        i++;
    }
    BasepairStateOP start (new BasepairState(start_conf->ends()[1]->state()));
    BasepairStateOP end (new BasepairState(start_conf->ends()[2]->state()));
    
    MotifTreeStateSearch search;
      MotifTreeStateSelectorOP selector ( new MotifTreeStateSelector(mts_libs, "round_robin"));
    //MotifTreeStateSelectorOP selector = custom_selector("fmn");
    search.base_beads(centers);
    search.set_numeric_option("max_size", 130);
    search.set_numeric_option("accept_ss_score", 0);
    search.set_numeric_option("max_n_solutions", 1);
    search.set_numeric_option("save_solutions", 1);
    search.set_numeric_option("accept_score", 10);

    MotifTreeStateSearchScorerOP scorer ( new MTSS_Astar() );
    MotifTreeStateSearchSolutionOPs solutions = search.search(start, end, selector, scorer);
    std::cout << solutions.size() << std::endl;

    
    i = 0;
    std::stringstream ss;
    for (auto const & s : solutions) {
        for (auto const & n : s->path() ) {
            std::cout << n->mts()->name() << " ";
        }
        std::cout << std::endl;
        ss << "solutions." << i << ".pdb";
        //s->to_mtst().to_pdb(ss.str());
        PoseOP sol = s->to_mtst().to_pose();
        mt.add_motif(sol);
        mt.write_pdbs();
        ss.str("");
        i++;
    }

    
    
    return 1;
}

int build_ribozyme_final_model() {
    
    Strings motifs = get_lines_from_file(base_dir() + "/rnamake/lib/RNAMake/motif_tree_path_refiner/hammerhead_motifs.str");
    ResidueTypeSet rts;
    MotifOP hammerhead_tc ( new Motif(motifs[0], rts));
    MotifOP hammerhead_core ( new Motif(motifs[1], rts));
    MotifOP fmn (new Motif(motifs[2], rts));
    
    MotifTree mt;
    MotifLibrary mlib ( HELIX );
    MotifTreeNodeOP node = mt.add_motif(hammerhead_core, NULL, 2);
    mt.add_motif(mlib.get_motif("HELIX.LE.20"));
    mt.add_motif(hammerhead_tc, NULL, 1);
    mt.add_motif(fmn, node);
    
    PoseOP start_conf = mt.to_pose();
    MotifTreeStateOP start_mts = motif_to_state(start_conf, 0);
    
    MotifTreeStateLibraryOPs mts_libs(2);
    mts_libs[0] = MotifTreeStateLibraryOP ( new MotifTreeStateLibrary(HELIX, 5));
    mts_libs[1] = MotifTreeStateLibraryOP ( new MotifTreeStateLibrary(TWOWAY));
    
    Strings names = split_str_by_delimiter("HELIX.LE.4-0-0-0-0-1-1 TWOWAY.2PXV.1-0-0-1-0-0-0 HELIX.LE.5-0-0-0-0-1-1 TWOWAY.1S72.49-0-0-0-0-1-1 HELIX.LE.8-0-0-0-0-1-1 TWOWAY.3P59.1-0-0-0-0-1-1 HELIX.LE.2-0-0-0-0-1-0 TWOWAY.2ZY6.0-0-0-0-0-1-1 HELIX.LE.4-0-0-0-0-1-1", " ");
    
    MotifTreeStateTree mtst;
    mtst.add_state(start_mts, NULL);
    mtst.sterics(0);
    
    int pos = 0;
    int i = -1;
    
    for(auto const & name : names) {
        i++;
        if ( i % 2 == 0)  { pos = 0; }
        else              { pos = 1; }
        
        if(i == 1) {
            mtst.add_state(mts_libs[pos]->get_state(name), NULL, 1);
        }
        else {
            mtst.add_state(mts_libs[pos]->get_state(name), NULL);
        }
        
    }

    mtst.to_pdb("test_2.pdb");
    
    return 1;
}


int main(int argc, const char * argv[]) {
    // insert code here...
    //test_path_refiner_ttr();
    //test_path_refiner_tecto(argc, argv);
    //test_converter();
    //test_sequence_designer();
    optimize(std::stoi(argv[1]));
    //build_ribozyme();
    //build_ribozyme_final_model();
    return 0;
}






































