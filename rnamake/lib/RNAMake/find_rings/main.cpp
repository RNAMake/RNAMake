//
//  main.cpp
//  find_rings
//
//  Created by Joseph Yesselman on 3/7/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include <map>
#include <sstream>
#include "types.h"
#include "chain.h"
#include "residue_type_set.h"
#include "settings.h"
#include "motif.h"
#include "motif_library.h"
#include "motif_tree_state.h"
#include "motif_tree_state_library.h"
#include "motif_tree_state_tree.h"
#include "cartesian_product.h"


float
inline
calculate_frame_score (
    BasepairStateOP const & target,
    BasepairStateOP const & current) {
    
    return current->d().distance(target->d()) + current->r().difference(target->r());
    
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

int test_cartesian_product() {
    Ints int_test(3);
    int_test[0] = 0; int_test[1] = 1; int_test[2] = 2;
    std::vector<Ints> values(3);
    values[0] = int_test;
    values[1] = int_test;
    values[2] = int_test;
    
    CartesianProduct<int> test_product(values);
    while(! test_product.end() ) {
        Ints current = test_product.next();
        for(auto const & c : current) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
    }
    
    return 1;
}

ChainOPs
get_chains_for_end(
    BasepairOP const & end,
    MotifOP const & m) {
    
    ChainOPs chains;
    for (auto const & c : m->chains()) {
        for (auto const & r : end->residues() ) {
            if(c->first() == r) { chains.push_back(c); }
            if(c->last()  == r) { chains.push_back(c); }
        }
    }
    
    return chains;
}


std::map<int, Ints>
get_different_chains(
    MotifOP const & m) {
    
    std::map<int, Ints> different_chains;
    int i = -1, j = -1;
    int same = 0;
    for (auto const & end1 : m->ends()) {
        i++;
        j = -1;
        ChainOPs end1_chains = get_chains_for_end(end1, m);
        for (auto const & end2 : m->ends()) {
            j++;
            if (i >= j) { continue; }
            same = 0;
            ChainOPs end2_chains = get_chains_for_end(end2, m);
            
            for (auto const & c1 : end1_chains) {
                for(auto const & c2 : end2_chains) {
                    if(c1 == c2) { same = 1; break; }
                }
            }
            
            if(same) { continue; }
            
            if(different_chains.find(i) == different_chains.end()) { different_chains[i] = Ints(); }
            if(different_chains.find(j) == different_chains.end()) { different_chains[j] = Ints(); }
            
            different_chains[i].push_back(j);
            different_chains[j].push_back(i);
            
        }
    }
    
    return different_chains;

}

MotifTreeStateOPs
get_tc_states(
    MotifOP const & m) {
 
    MotifTreeStateOPs all_tc_states;
    for (int i = 0; i < m->ends().size(); i++) {
        for(int j = 0; j < 2; j++) {
            MotifTreeStateOP state = motif_to_state(m, i, j);
            all_tc_states.push_back(state);
        }
    }
    
    return all_tc_states;
}


void
find_rings(
    MotifOP const & m1,
    MotifOP const & m2,
    int i,
    int j,
    MotifTreeStateLibrary const & mts_lib) {
    
    std::map<int, Ints> different_chains1 = get_different_chains(m1);
    std::map<int, Ints> different_chains2 = get_different_chains(m2);
    MotifTreeStateOPs tc_states_1 = get_tc_states(m1);
    MotifTreeStateOPs tc_states_2 = get_tc_states(m2);
    MotifTreeStateOPs helices;
    for (auto const & mts : mts_lib.motif_tree_states()) {
        if (mts->size() < 10) { continue; }
        if (mts->size() > 30) { continue; }
        helices.push_back(mts);
    }
    
    std::vector<MotifTreeStateOPs> topology(4);
    topology[0] = tc_states_1;
    topology[1] = helices;
    topology[2] = tc_states_2;
    topology[3] = helices;
    
    CartesianProduct<MotifTreeStateOP> topology_iterator(topology);
    MotifTreeStateOPs c;
    
    MotifTreeStateTree mtst;
    MotifTreeStateNodeOP node;
    BasepairStateOP final_end, current;
    Ints end_indexes, end_indexes2;
    float dist, dist2;
    int count = 0;
    std::stringstream ss;
    
    while(!topology_iterator.end()) {
        c = topology_iterator.next();
        mtst.add_state(c[0], NULL);
        end_indexes = different_chains1 [ c[0]->start_index() ];
        for( auto const & ei : end_indexes ) {
            node = mtst.add_state(c[1], NULL, ei);
            if(node == NULL) { continue; }
            
            node = mtst.add_state(c[2], NULL);
            
            if(node == NULL) {
                mtst.remove_node();
                continue;
            }
            
            final_end = mtst.nodes()[0]->active_states()[0];
            end_indexes2 = different_chains2 [ c[2]->start_index() ];
            for (auto const & ei2 : end_indexes2) {
                node = mtst.add_state(c[3], NULL, ei2);
                if(node == NULL) { continue; }
                current = mtst.last_node()->active_states()[0];
                dist = calculate_frame_score(final_end, current);
                current->flip();
                dist2 = calculate_frame_score(final_end, current);
                current->flip();
                if(dist2 < dist) { dist = dist2; }
                if( dist < 10) {
                    ss << "ring." << i << "." << j << "." << count << ".pdb";
                    count ++;
                    mtst.to_pdb(ss.str());
                    ss.str("");
                    return;
                }
                mtst.remove_node();
            }
            mtst.remove_node();
            mtst.remove_node();
        }
        
        mtst.remove_node();
    }

    

    
}



int main(int argc, const char * argv[]) {
    // insert code here...
    
    Strings lines = get_lines_from_file(base_dir() + "/rnamake/lib/RNAMake/simulate_tectos/tetraloop.str");
    ResidueTypeSet rts;
    MotifOP gaaa_motif ( new Motif(lines[1], rts));
    MotifLibrary mlib (TCONTACT);
    mlib.load_all(87);
    MotifOPs motifs = mlib.motifs();
    //motifs.push_back(gaaa_motif);
    
    MotifTreeStateLibrary mts_lib (HELIX);
    int i = -1, j = -1;
    
    for(auto const & m1 : motifs) {
        i++;
        std::cout << i << std::endl;
        find_rings(gaaa_motif, m1, 0, i+0, mts_lib);
    }
    
    /*for(auto const & m1 : motifs) {
        j = -1;
        i++;
        if(i < 74) { continue;}
        
        for(auto const & m2 : motifs) {
            j++;
            if(i < 75 && j < 87 ) { continue; }
            if( i > j) { continue; }
            std::cout << i << " " << j << std::endl;
            find_rings(m1, m2, i, j, mts_lib);
        }
    }*/
    
    
    exit(0);
    
    return 0;
}






















