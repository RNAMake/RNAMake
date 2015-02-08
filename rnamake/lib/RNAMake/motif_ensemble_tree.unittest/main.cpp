//
//  main.cpp
//  motif_ensemble_tree.unittest
//
//  Created by Joseph Yesselman on 2/6/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include <time.h>
#include <random>

#include "motif_ensemble.h"
#include "motif_ensemble_tree.h"
#include "motif_tree_state.h"
#include "motif_tree_state_tree.h"

inline float GenXORShift32(void)
{
    static unsigned seed = 2463534242U;
    
    seed ^= (seed << 5);
    seed ^= (seed >> 13);
    seed ^= (seed << 6);
    
    return seed * (1.0f / 4294967295.0f);
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
    Strings steps = split_str_by_delimiter("GC=GC,GC=UA,UA=AU,AU=AU,AU=GC,GC=GU,GU=GC,GC=CG,CG=AU,AU=GC,GC=AU,AU=GC,GC=GC", ",");
    //Strings steps = split_str_by_delimiter("GC=GC,GC=UA,UA=AU,AU=AU,AU=GC,GC=GC,GC=GC,GC=CG,CG=AU,AU=GC,GC=AU,AU=GC,GC=GC", ",");
    for(auto const & step : steps) {
        MotifEnsemble me( step, 0, 0);
        met.add_ensemble(me);
    }
    
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<float> rdist(0,1);

    
    MotifTreeStateTree mtst = met.get_mtst();
    MotifEnsembleTreeNodeOP met_node;
    MotifTreeStateNodeOP mts_node;
    MotifState ms = met.nodes()[0]->motif_ensemble().motif_states()[0];
    float kB =  1.3806488e-1 ;
    float kBT =  kB * 298.15 ;
    float score;
    int node_num = 0;
    int i = 0, nsteps =10000000, result = 0;
    float cpop;
    float dist;
    std::ofstream out;
    std::vector<float> bins(120);
    float bin_size = 0.25;
    float min_dist = 20;
    float max_dist = 50;
    float bin_pos;
    int pos;
    srand(time(NULL));
    while (i < nsteps) {
        if( i % 10000 == 0 ) { std::cout << i << std::endl; }
        node_num = 1 + (mtst.nodes().size()-1)*rdist(mt);
        if(node_num == 0 ) { continue; }
        met_node = met.nodes()[ node_num ];
        mts_node = mtst.nodes()[ node_num ];
        cpop = met_node->motif_ensemble().get_state(mts_node->mts()->name()).population;
        pos = (int)((met_node->motif_ensemble().motif_states().size()-1)*rdist(mt));
        ms = met_node->motif_ensemble().get_state(pos);
        if ( ms.population < cpop) {
            result = mtst.replace_state(mts_node, ms.mts);
            dist =  mtst.nodes()[0]->states()[0]->sugars()[0].distance(mtst.nodes().back()->states()[1]->sugars()[1]);
            bin_pos = (int)((dist - min_dist) / bin_size);
            bins[bin_pos] += 1;
            i++;
            
        }
        
        score = exp((cpop - ms.population)/kBT);
        if( rdist(mt) < score) {
            result = mtst.replace_state(mts_node, ms.mts);
            dist =  mtst.nodes()[0]->states()[0]->sugars()[0].distance(mtst.nodes().back()->states()[1]->sugars()[1]);
            bin_pos = (int)((dist - min_dist) / bin_size);
            bins[bin_pos] += 1;
            i++;
        }
    }
    
    //MotifTree mt = mtst.to_motiftree();
    i = 0;
    for (auto const & b : bins) {
        std::cout << min_dist + i*bin_size << " " << b << std::endl;
        i++;
    }
    
    
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
    //Strings chip_steps = split_str_by_delimiter("GC=AU AU=AU AU=GC GC=AU AU=UA UA=CG CG=CG CG=UA UA=CG CG=GC", " ");
    //Strings flow_steps = split_str_by_delimiter("GC=AU AU=AU AU=GC GC=UA UA=AU AU=CG CG=CG CG=GC GC=GU GU=GC", " ");
    Strings chip_steps = split_str_by_delimiter("GC=AU AU=AU AU=GC GC=AU AU=UA UA=CG CG=CG CG=UG UG=GU GU=GC", " ");

    MotifEnsembleTree met;
    met.add_ensemble(MotifEnsemble("AU=GC", 0, 0));
    met.add_ensemble(mts_to_me(ggaa_state));
    met.add_ensemble(MotifEnsemble(flow_steps[0], 0, 1), NULL, 2);
    for (int i = 1; i < flow_steps.size(); i++) {
        met.add_ensemble(MotifEnsemble(flow_steps[i], 0, 1));
    }
    met.add_ensemble(mts_to_me(gaaa_state));
    met.add_ensemble(MotifEnsemble(chip_steps[0], 0, 1), NULL, 1);
    for (int i = 1; i < chip_steps.size(); i++) {
        met.add_ensemble(MotifEnsemble(chip_steps[i], 0, 1));
    }

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<float> dist(0,1);

    MotifTreeStateTree mtst = met.get_mtst();
    MotifEnsembleTreeNodeOP met_node;
    MotifTreeStateNodeOP mts_node;
    MotifState ms = met.nodes()[0]->motif_ensemble().motif_states()[0];
    float kB =  1.3806488e-1 ;
    float kBT =  kB * 298.15 ;
    float score;
    int node_num = 0;
    int i = 0, nsteps =10000000, result = 0;
    float cpop;
    float frame_score;
    BasepairStateOP target = mtst.nodes()[2]->states()[0];
    BasepairStateOP target_flip (new BasepairState(target->copy()));
    target_flip->flip();
    int under_cutoff = 0;
    float cutoff = 5.0f;
    float diceroll = 0.0f;
    int pos = 0;
    
    srand(unsigned(time(NULL)));
    while (i < nsteps) {
        //if( i % 10000 == 0 ) { std::cout << under_cutoff <<" " << i << std::endl; }
        node_num = (int)((met.nodes().size()-1)*dist(mt));
        if(node_num == 0) { continue; }
        //node_num = 1 + (int)((met.nodes().size()-2)*GenXORShift32());
        met_node = met.nodes()[ node_num ];
        mts_node = mtst.nodes()[ node_num ];
        cpop = met_node->motif_ensemble().get_state(mts_node->mts()->name()).population;
        pos = (int)((met_node->motif_ensemble().motif_states().size()-1)*dist(mt));
        ms = met_node->motif_ensemble().get_state(pos);
        if ( ms.population < cpop) {
            result = mtst.replace_state(mts_node, ms.mts);
            frame_score = frame_distance(mtst.last_node()->states()[1], target, target_flip);
            if(frame_score < cutoff) { under_cutoff++; }
            i++;
            continue;

        }
        
        score = exp((cpop - ms.population)/kBT);
        diceroll = dist(mt);
        if( diceroll < score) {
            result = mtst.replace_state(mts_node, ms.mts);
            frame_score = frame_distance(mtst.last_node()->states()[1], target, target_flip);
            if(frame_score < cutoff) { under_cutoff++; }
            i++;
        }
    }

    std::cout << under_cutoff << " " << nsteps << std::endl;
    
    
    //MotifTreeStateTree mtst = met.get_mtst();
    //MotifTree mt = mtst.to_motiftree();
    //mt.write_pdbs();
    
    return 1;
}

int main(int argc, const char * argv[]) {
    //if (test_creation_me() == 0)          { std::cout << "test_creation_me failed" << std::endl;  }
    //if (test_creation_me_tree() == 0)     { std::cout << "test_creation_me_tree failed" << std::endl;  }
    //if (test_add_ensemble() == 0)         { std::cout << "test_add_ensemble failed" << std::endl;  }
    //if (test_sample() == 0)               { std::cout << "test_sample failed" << std::endl;  }
    String flow_string = "GC=AU,AU=AU,AU=GC,GC=UA,UA=AU,AU=CG,CG=CG,CG=GC,GC=AU,AU=GC";
    String chip_string = "GC=AU,AU=AU,AU=GC,GC=AU,AU=UA,UA=CG,CG=CG,CG=UG,UG=GU,GU=GC";
    
    
    
    test_tecto();
    
    return 0;
}


















