//
//  sample_helix.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 6/7/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "sample_helix/sample_helix.hpp"
#include "base/backtrace.hpp"
#include "secondary_structure/secondary_structure_parser.h"
#include "motif_data_structures/motif_tree.h"
#include "motif_data_structures/motif_state_ensemble_tree.h"
#include "resources/resource_manager.h"
#include "thermo_fluctuation/thermo_fluc_sampler.h"

SampleHelixApp::SampleHelixApp() : Application(),
        sampler_(ThermoFlucSampler())
{}


// application setups functions ////////////////////////////////////////////////////////////////////

void
SampleHelixApp::setup_options() {
    add_option("seq", "", OptionType::STRING, true);
    add_option("s", 1000000, OptionType::INT);
    add_option("output_freq", 100, OptionType::INT);
    add_option("output_filename", "test.out", OptionType::STRING);


}

void
SampleHelixApp::parse_command_line(
        int argc,
        const char ** argv) {
    Application::parse_command_line(argc, argv);
}

void
SampleHelixApp::run() {
    auto sequence = get_string_option("seq");
    auto structure = _generate_structure(sequence);

    auto motifs = get_motifs_from_seq_and_ss(sequence, structure);
    auto mt = std::make_shared<MotifTree>();
    for(auto const & m : motifs) {
        mt->add_motif(m);
    }

    auto mset = std::make_shared<MotifStateEnsembleTree>(mt);

    sampler_.setup(mset);

    std::ofstream out;
    out.open(get_string_option("output_filename"));

    int output_freq = get_int_option("output_freq");
    for(int i = 0; i < get_int_option("s"); i++) {
        sampler_.next();
        if(i % output_freq == 0) {
            out << sampler_.mst()->last_node()->data()->get_end_state(1)->d().to_str() << "&";
            out << sampler_.mst()->last_node()->data()->get_end_state(1)->r().to_str() << "&";
            out << sampler_.mst()->topology_to_str() << std::endl;
        }
    }
    out.close();


}

String
SampleHelixApp::_generate_structure(
        String const & seq) {
    auto structure = String("");
    auto spl = split_str_by_delimiter(seq, "&");

    if(spl.size() != 2) {
        throw std::runtime_error(
                "sequence must be composed of two strands");
    }

    for(int i = 0; i < spl[0].size(); i++) { structure += "("; }
    structure += "&";
    for(int i = 0; i < spl[0].size(); i++) { structure += ")"; }
    return structure;
}

MotifOPs
SampleHelixApp::get_motifs_from_seq_and_ss(
        String const & seq,
        String const & ss) {

    auto parser = sstruct::SecondaryStructureParser();
    auto ss_motifs = parser.parse_to_motifs(seq, ss);
    auto motifs = MotifOPs();

    auto start = 0;
    auto motif = MotifOP(nullptr);
    for(auto const & m : ss_motifs) {
        //basepair step
        if(m->mtype() == MotifType::HELIX) {
            motif = RM::instance().bp_step(m->end_ids()[0]);
            motifs.push_back(motif);
        }
        else {
            throw std::runtime_error("only helices are allowed");
        }
    }

    return motifs;
}


int main(int argc, const char * argv[]) {
    //must add this for all apps!
    std::set_terminate(print_backtrace);

    auto app = SampleHelixApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    /*auto seq = String("GGCCCUCAAGGG&CCCUUGAGGGCC");
    auto ss  = String("((((((((((((&))))))))))))");
    
    sstruct::SS_Tree ss_tree(seq, ss);
    sstruct::SS_TreeNodeOP last_node = nullptr;

    auto mt = std::make_shared<MotifTree>();
    
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
    
    auto mf = MotifFactory();
    auto m = mf.motif_from_file("/Users/josephyesselman/Downloads/sele.pdb");
    m->block_end_add(-1);
    m->ends()[0]->flip();
    m->mtype(HAIRPIN);
    ResourceManager::getInstance().register_motif(m);
    mt->add_motif(m);

    auto mset = std::make_shared<MotifStateEnsembleTree>();
    mset->setup_from_mt(mt);
    
    auto sampler = ThermoFlucSampler();
    sampler.setup(mset);
    sampler.to_pdb("start.pdb");
    
    for(int i = 0; i < 10000000; i++) {
        try {
            sampler.next();
        } catch(...) { }
            
        if(i % 10000 == 0 ) {
            try {
                sampler.to_pdb("sampled."+std::to_string(i)+".pdb");
            } catch(...) { }
        }
        
    }*/

    
    return 0;
}