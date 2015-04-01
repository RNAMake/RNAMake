//
//  main.cpp
//  motif_tree_state.unittest
//
//  Created by Joseph Yesselman on 2/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "motif_tree_state.h"
#include "motif_tree_state_library.h"
#include "motif_tree_state_node.h"
#include "motif_tree_state_node_aligner.h"
#include "motif_tree_state_tree.h"

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

int
test_creation() {
    MotifTreeStateLibrary mts_lib ( TWOWAY );
    return 1;
}

int
test_creation_node() {
    MotifTreeStateLibrary mts_lib ( TWOWAY );
    MotifTreeStateNode mtn ( mts_lib.motif_tree_states()[0], 0, NULL, 0);
    Ints avail_ends = mtn.available_ends();
    return 1;
}

int
test_add_child() {
    MotifTreeStateLibrary mts_lib ( TWOWAY );
    MotifTreeStateNodeOP mtn1 (new MotifTreeStateNode( mts_lib.motif_tree_states()[0], 0, NULL, 0));
    MotifTreeStateNodeOP mtn2 (new MotifTreeStateNode( mts_lib.motif_tree_states()[0], 0, mtn1, 0));
    Ints indices = mtn1->available_ends();
    mtn1->add_child(mtn2, indices[0]);
    indices = mtn1->available_ends();
    int parent_end_index = mtn2->parent_end_index();
    BasepairStateOP parent_end = mtn2->parent_end();
    return 1;
}

int
test_replace_mts() {
    MotifTreeStateLibrary mts_lib ( TWOWAY );
    MotifTreeStateNodeOP mtn1 (new MotifTreeStateNode( mts_lib.motif_tree_states()[0], 0, NULL, 0));
    mtn1->replace_mts(mts_lib.motif_tree_states()[1]);
    return 1;
}

int
test_creation_mtst() {
    MotifTreeStateTree mtst;
    return 1;
}

int
test_add_state() {
    Strings lines = get_lines_from_file("test_add_state.dat");
    MotifTreeStateLibrary mts_lib ( TWOWAY );
    for(auto const & l : lines) {
        Strings spl = split_str_by_delimiter(l, " ");
        MotifTreeStateTree mtst;
        for (int i = 1; i < spl.size(); i++) {
            MotifTreeStateOP mts = mts_lib.get_state(spl[i]);
            mtst.add_state(mts, NULL, NULL);
        }
        
        if (spl.size() != mtst.nodes().size()) {
            std::cout << spl.size() << " " << mtst.nodes().size() << " " << std::endl;
        }
    }
    return 1;
}

int
test_compare_last_node() {
    Strings lines = get_lines_from_file("test_compare_last_node.dat");
    MotifTreeStateLibrary mts_lib ( TWOWAY );
    for(auto const & l : lines) {
        Strings spl = split_str_by_delimiter(l, "|");
        Strings names = split_str_by_delimiter(spl[0], " ");
        BasepairState target = str_to_basepairstate(spl[1]);
        MotifTreeStateTree mtst;
        for (int i = 1; i < names.size(); i++) {
            MotifTreeStateOP mts = mts_lib.get_state(names[i]);
            mtst.add_state(mts, NULL, NULL);
        }
        int last_index = mtst.nodes().back()->available_ends()[0];
        BasepairStateOP state = mtst.nodes().back()->states()[last_index];
        float dist = state->d().distance(target.d());
        if(dist > 0.1) {
            std::cout << dist << std::endl;
        }
        float matrix_dist = state->r().difference(target.r());
        if (matrix_dist > 0.1) {
            std::cout << matrix_dist << std::endl;
        }
    }

    return 1;
}

int
test_to_motiftree() {
    Strings lines = get_lines_from_file("test_add_state.dat");
    MotifTreeStateLibrary mts_lib ( TWOWAY );
    for(int k = 0; k < 10000; k++) {
    for(auto const & l : lines) {
        Strings spl = split_str_by_delimiter(l, " ");
        MotifTreeStateTree mtst;
        for (int i = 1; i < spl.size(); i++) {
            MotifTreeStateOP mts = mts_lib.get_state(spl[i]);
            mtst.add_state(mts, NULL, NULL);
        }
        MotifTree mt = mtst.to_motiftree();
        std::cout << k << std::endl;
        //mt.write_pdbs();
        //break;

    }
    }
    return 1;
}

int
test_replace_state() {
    Strings lines = get_lines_from_file("test_add_state.dat");
    MotifTreeStateLibrary mts_lib ( TWOWAY );
    Strings spl = split_str_by_delimiter(lines[0], " ");
    MotifTreeStateTree mtst;
    for (int i = 1; i < spl.size(); i++) {
        MotifTreeStateOP mts = mts_lib.get_state(spl[i]);
        mtst.add_state(mts, NULL, NULL);
    }
    MotifTree mt = mtst.to_motiftree();
    //mt.write_pdbs("org");
    int result = 0;
    for (auto const & mts : mts_lib.motif_tree_states()) {
        result = mtst.replace_state(mtst.nodes()[10], mts);
        if(result == 1) { break; }
    }
    mt = mtst.to_motiftree();
    //mt.write_pdbs();

    return 1;
}



int main(int argc, const char * argv[]) {
    //if (test_creation() == 0)          { std::cout << "test_creation failed" << std::endl;  }
    //if (test_creation_node() == 0)     { std::cout << "test_creation_node failed" << std::endl;  }
    //if (test_add_child() == 0)         { std::cout << "test_add_child failed" << std::endl;  }
    //if (test_replace_mts() == 0)       { std::cout << "test_replace_mts failed" << std::endl;  }
    //if (test_creation_mtst() == 0)     { std::cout << "test_creation_mtst failed" << std::endl;  }
    //if (test_add_state() == 0)        { std::cout << "test_add_state failed" << std::endl;  }
    //if (test_compare_last_node() == 0) { std::cout << "test_compare_last_node failed" << std::endl;  }
    if (test_to_motiftree() == 0)      { std::cout << "test_to_motiftree failed" << std::endl;  }
    //if (test_replace_state() == 0)     { std::cout << "test_replace_state failed" << std::endl;  }

    return 0;
}










