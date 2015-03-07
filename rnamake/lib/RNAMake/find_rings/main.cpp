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

int main(int argc, const char * argv[]) {
    // insert code here...
    
    Strings lines = get_lines_from_file(base_dir() + "/rnamake/lib/RNAMake/simulate_tectos/tetraloop.str");
    ResidueTypeSet rts;
    MotifOP gaaa_motif ( new Motif(lines[1], rts));
    
    std::map<int, Ints> different_chains;
    int i = -1, j = -1;
    int same = 0;
    for (auto const & end1 : gaaa_motif->ends()) {
        i++;
        j = -1;
        ChainOPs end1_chains = get_chains_for_end(end1, gaaa_motif);
        for (auto const & end2 : gaaa_motif->ends()) {
            j++;
            if (i >= j) { continue; }
            same = 0;
            ChainOPs end2_chains = get_chains_for_end(end2, gaaa_motif);
            
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
    
    MotifTreeStateOPs all_tc_states;
    for (i = 0; i < gaaa_motif->ends().size(); i++) {
        for(j = 0; j < 2; j++) {
            MotifTreeStateOP gaaa_state = motif_to_state(gaaa_motif, i, j);
            all_tc_states.push_back(gaaa_state);
        }
    }
    
    MotifTreeStateLibrary mts_lib (HELIX);
    
    std::vector<MotifTreeStateOPs> topology(4);
    topology[0] = all_tc_states;
    topology[1] = mts_lib.motif_tree_states();
    topology[2] = all_tc_states;
    topology[3] = mts_lib.motif_tree_states();
    
    CartesianProduct<MotifTreeStateOP> topology_iterator(topology);
    MotifTreeStateOPs c;
    
    MotifTreeStateTree mtst;
    MotifTreeStateNodeOP node;
    BasepairStateOP final_end, current;
    Ints end_indexes, end_indexes2;
    float dist, dist2;
    int count = -1;
    std::stringstream ss;
    while(!topology_iterator.end()) {
        c = topology_iterator.next();
        count++;
        if(count % 100 == 0) { std::cout << count << std::endl; }
        mtst.add_state(c[0], NULL);
        end_indexes = different_chains [ c[0]->start_index() ];
        for( auto const & ei : end_indexes ) {
            node = mtst.add_state(c[1], NULL, ei);
            if(node == NULL) { continue; }
            
            node = mtst.add_state(c[2], NULL);
            
            if(node == NULL) {
                mtst.remove_node();
                continue;
            }
            
            final_end = mtst.nodes()[0]->active_states()[0];
            end_indexes2 = different_chains [ c[2]->start_index() ];
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
                    std::cout << dist << std::endl;
                    ss << "success." << count << ".pdb";
                    try {
                        mtst.to_pdb(ss.str());
                    }
                    catch(String e) { std::cout << "caught" << std::endl;}
                    ss.str("");
                }
                mtst.remove_node();
            }
            mtst.remove_node();
            mtst.remove_node();
        }
        
        mtst.remove_node();
    }
    
    
    return 0;
}






















