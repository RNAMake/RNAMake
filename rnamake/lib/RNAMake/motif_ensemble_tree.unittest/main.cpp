//
//  main.cpp
//  motif_ensemble_tree.unittest
//
//  Created by Joseph Yesselman on 2/6/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include <time.h>
#include "motif_ensemble.h"
#include "motif_ensemble_tree.h"
#include "motif_tree_state.h"
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
test_creation_me() {
    MotifEnsemble me ("GC=GC", 0, 0);
    return 1;
}

int
test_creation_me_tree() {
    MotifEnsemble me ("GC=GC", 0, 0);
    MotifEnsembleTree met;
    MotifEnsembleTree met2 (me);
    return 1;
}

int
test_add_ensemble() {
    MotifEnsembleTree met;
    MotifEnsemble me ("GC=GC", 0, 0);
    MotifEnsemble me2 ("GC=GC", 0, 0);
    met.add_ensemble(me);
    met.add_ensemble(me2);
    MotifTreeStateTree mtst = met.get_mtst();
    for(auto const & n : mtst.nodes()) {
     //   std::cout << n->mts()->name() << std::endl;
    }
    MotifTree mt = mtst.to_motiftree();
    return 1;
}

int
test_sample() {
    MotifEnsembleTree met;
    Strings steps = split_str_by_delimiter("GC=GC,GC=UA,UA=AU,AU=AU,AU=GC,GC=GC,GC=GU,GU=CG,CG=AU,AU=GC,GC=AU,AU=GC,GC=GC", ",");
    for(auto const & step : steps) {
        MotifEnsemble me( step, 0, 0);
        met.add_ensemble(me);
    }
    MotifTreeStateTree mtst = met.get_mtst();
    MotifEnsembleTreeNodeOP met_node;
    MotifTreeStateNodeOP mts_node;
    MotifState ms = met.nodes()[0]->motif_ensemble().motif_states()[0];
    float kB =  1.3806488e-1 ;
    float kBT =  kB * 298.15 ;
    float score;
    int node_num = 0;
    int i = 0, nsteps =100000, result = 0;
    float cpop;
    float dist;
    std::ofstream out;
    out.open("test_gu.dat");
    srand(time(NULL));
    while (i < nsteps) {
        if( i % 10000 == 0 ) { std::cout << i << std::endl; }
        node_num = 1 + rand() % mtst.nodes().size()-1;
        if(node_num == 0 ) { continue; }
        met_node = met.nodes()[ node_num ];
        mts_node = mtst.nodes()[ node_num ];
        cpop = met_node->motif_ensemble().get_state(mts_node->mts()->name()).population;
        ms = met_node->motif_ensemble().get_random_state();
        if ( ms.population < cpop) {
            result = mtst.replace_state(mts_node, ms.mts);
            if( i % 100 == 0) {
                dist =  mtst.nodes()[0]->states()[0]->sugars()[1].distance(mtst.nodes().back()->states()[1]->sugars()[1]);
                out << dist << std::endl;
            }
        }
        
        score = exp((cpop - ms.population)/kBT);
        if( rand() % 1 < score) {
            result = mtst.replace_state(mts_node, ms.mts);
            if (i % 100) {
                dist =  mtst.nodes()[0]->states()[0]->sugars()[1].distance(mtst.nodes().back()->states()[1]->sugars()[1]);
                out << dist << std::endl;
            }
        }
        i++;
    }
    out.close();
    
    //MotifTree mt = mtst.to_motiftree();
    
    
    return 1;
}

MotifEnsemble
mts_to_me(MotifTreeStateOP const & mts) {
    MotifEnsemble me;
    MotifState ms (mts, 1.0);
    me.add_motif_state(ms);
    return me;
}

int
test_tecto() {
    Strings lines = get_lines_from_file("tetraloop.str");
    ResidueTypeSet rts;
    MotifOP ggaa_motif ( new Motif(lines[0], rts));
    MotifOP gaaa_motif ( new Motif(lines[1], rts));
    MotifTreeStateOP ggaa_state = motif_to_state(ggaa_motif, 1, 1);
    MotifTreeStateOP gaaa_state = motif_to_state(gaaa_motif, 0, 1);
    Strings flow_steps = split_str_by_delimiter("GC=AU AU=AU AU=GC GC=UA UA=AU AU=CG CG=CG CG=GC GC=AU AU=GC", " ");
    Strings chip_steps = split_str_by_delimiter("GC=AU AU=AU AU=GC GC=AU AU=UA UA=CG CG=CG CG=UA UA=CG CG=GC", " ");
    MotifEnsembleTree met;
    met.add_ensemble(MotifEnsemble("AU=GC", 0, 0));
    met.add_ensemble(mts_to_me(ggaa_state));
    met.add_ensemble(MotifEnsemble(flow_steps[0], 0, 0), NULL, 2);
    for (int i = 1; i < flow_steps.size(); i++) {
        met.add_ensemble(MotifEnsemble(flow_steps[i], 0, 0));
    }
    met.add_ensemble(mts_to_me(gaaa_state));
    met.add_ensemble(MotifEnsemble(chip_steps[0], 0, 0), NULL, 1);
    for (int i = 1; i < chip_steps.size(); i++) {
        met.add_ensemble(MotifEnsemble(chip_steps[i], 0, 0));
    }
    
    MotifTreeStateTree mtst = met.get_mtst();
    MotifTree mt = mtst.to_motiftree();
    mt.write_pdbs();
    
    return 1;
}

int main(int argc, const char * argv[]) {
    //if (test_creation_me() == 0)          { std::cout << "test_creation_me failed" << std::endl;  }
    //if (test_creation_me_tree() == 0)     { std::cout << "test_creation_me_tree failed" << std::endl;  }
    //if (test_add_ensemble() == 0)         { std::cout << "test_add_ensemble failed" << std::endl;  }
    //if (test_sample() == 0)               { std::cout << "test_sample failed" << std::endl;  }
    test_tecto();
    
    return 0;
}


















