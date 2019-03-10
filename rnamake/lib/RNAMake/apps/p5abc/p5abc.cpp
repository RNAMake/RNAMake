//
//  p5abc.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 6/23/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "p5abc.hpp"

#include "base/cl_option.h"
#include "resources/resource_manager.h"
#include "motif/motif_factory.h"
#include "thermo_fluctuation/thermo_fluc_sampler.h"
#include "thermo_fluctuation/thermo_fluc_scorer.h"


Options
parse_command_line(
    int argc,
    const char ** argv) {
    
    CL_Options cl_opts;
    cl_opts.add_option("seq", "", STRING_TYPE,
                       "UAGAGUUUC&GGGACUCUG", false);
    cl_opts.add_option("s", "steps", FLOAT_TYPE, "1000000", false);
    
    return cl_opts.parse_command_line(argc, argv);
    
}


void
setup_motif() {
    auto base_path = base::base_dir() + "/rnamake/lib/RNAMake/apps/p5abc/";
    
    auto mf = MotifFactory();
    auto m = mf.motif_from_file(base_path+"p4p6_start.pdb");
    m->block_end_add(-1);
    m->ends()[0]->flip();
    
    auto m_clash_section = mf.motif_from_file(base_path+"p4p6_clash_section.pdb");

}


int main(int argc, const char * argv[]) {
    auto opts = parse_command_line(argc, argv);
    auto rts = ResidueTypeSet();

    
    auto seq = opts.option<String>("seq");
    auto ss  = String("(((((((((&)))))))))");
    
    
    //need to build from 3' to 5' based on where I am building from so much reverse teh sequence
    auto spl = split_str_by_delimiter(seq, "&");
    std::reverse(spl[0].begin(), spl[0].end());
    std::reverse(spl[1].begin(), spl[1].end());
    auto flip_sequence = spl[1] + "&" + spl[0];
    

    auto base_path = base::base_dir() + "/rnamake/lib/RNAMake/apps/p5abc/";
    auto m = file_to_motif(base_path+"start.motif");
    auto m_clash_section = file_to_motif(base_path+"clash_section.motif");
    
    ResourceManager::getInstance().register_motif(m);
    
    auto mt = std::make_shared<MotifTree>();
    mt->add_motif(m);
    
    sstruct::SS_Tree ss_tree(flip_sequence, ss);
    sstruct::SS_TreeNodeOP last_node = nullptr;
    
    int i = -1;
    for(auto const & n : ss_tree) {
        i++;
        if(i == 0) { continue; }
        if(n->data()->sequence() == "&&&") { continue; }
        if(last_node == nullptr ) { last_node = n; continue; }
        
        auto seq1 = last_node->data()->sequence();
        auto seq2 = n->data()->sequence();
        auto motif_name = String("");
        motif_name.push_back(seq1[0]); motif_name.push_back(seq1[2]);
        motif_name.push_back('=');
        motif_name.push_back(seq2[0]); motif_name.push_back(seq2[2]);
        
        last_node = n;
        mt->add_motif(ResourceManager::getInstance().get_motif(motif_name));
        
    }
    
    
    auto scorer = FrameScorer();
    auto ref_bp_end = m_clash_section->ends()[0]->state();
    
    auto mset = std::make_shared<MotifStateEnsembleTree>();
    mset->setup_from_mt(mt);
    
    auto sampler = ThermoFlucSampler();
    sampler.setup(mset);

    auto cut_off = 4.5f;
    auto bp_state = BasepairStateOP();
    auto score = 0.0f;
    auto count = 0;
    for(int i = 0; i < 1000000; i++) {
        try {
            sampler.next();
        } catch(...) { continue; }
        
        bp_state = sampler.mst()->last_node()->data()->cur_state->end_states()[1];
        score = scorer.score(bp_state, ref_bp_end);
        if(score < cut_off) {
            count += 1;
        }
        
    }
    
    std::cout << count << std::endl;
    
    return 0;
}
