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
#include "motif_data_structure/motif_tree.h"
#include "motif_data_structure/motif_state_ensemble_tree.h"
#include "resources/resource_manager.h"
#include "thermo_fluctuation/thermo_fluc_sampler.h"

SampleHelixApp::SampleHelixApp() : base::Application(),
        sampler_(thermo_fluctuation::ThermoFlucSampler())
{}


// application setups functions ////////////////////////////////////////////////////////////////////

void
SampleHelixApp::setup_options() {
    add_option("seq", "", base::OptionType::STRING, true);
    add_option("s", 1000000, base::OptionType::INT);
    add_option("output_freq", 100, base::OptionType::INT);
    add_option("output_filename", "test.out", base::OptionType::STRING);
    add_option("pdbs", false, base::OptionType::BOOL);



}

void
SampleHelixApp::parse_command_line(
        int argc,
        const char ** argv) {
    base::Application::parse_command_line(argc, argv);
}

void
SampleHelixApp::run() {
    auto output_pdbs = get_bool_option("pdbs");
    auto pdb_count = 0;
    auto sequence = get_string_option("seq");
    auto structure = _generate_structure(sequence);

    auto motifs = get_motifs_from_seq_and_ss(sequence, structure);
    auto mt = std::make_shared<motif_data_structure::MotifTree>();
    for(auto const & m : motifs) {
        mt->add_motif(m);
    }

    auto mset = std::make_shared<motif_data_structure::MotifStateEnsembleTree>(mt);

    sampler_.setup(mset);

    std::ofstream out;
    if(! output_pdbs) {
        out.open(get_string_option("output_filename"));
    }

    int output_freq = get_int_option("output_freq");
    for(int i = 0; i < get_int_option("s"); i++) {
        sampler_.next();
        if(i % output_freq == 0) {
            if(! output_pdbs) {
                out << sampler_.mst()->last_node()->data()->get_end_state(1)->d().to_str() << "&";
                out << sampler_.mst()->last_node()->data()->get_end_state(1)->r().to_str() << "&";
                out << sampler_.mst()->topology_to_str() << std::endl;
            }
            else {
                try {
                    sampler_.mst()->to_motif_tree()->to_pdb("test."+std::to_string(pdb_count)+".pdb", 1, 1);
                    pdb_count += 1;
                }
                catch(...) { }
            }
        }
    }
    out.close();


}

String
SampleHelixApp::_generate_structure(
        String const & seq) {
    auto structure = String("");
    auto spl = base::split_str_by_delimiter(seq, "&");

    if(spl.size() != 2) {
        throw std::runtime_error(
                "sequence must be composed of two strands");
    }

    for(int i = 0; i < spl[0].size(); i++) { structure += "("; }
    structure += "&";
    for(int i = 0; i < spl[0].size(); i++) { structure += ")"; }
    return structure;
}

motif::MotifOPs
SampleHelixApp::get_motifs_from_seq_and_ss(
        String const & seq,
        String const & ss) {

    auto parser = secondary_structure::Parser();
    auto ss_motifs = parser.parse_to_motifs(seq, ss);
    auto motifs = motif::MotifOPs();

    auto start = 0;
    auto motif = motif::MotifOP(nullptr);
    for(auto const & m : ss_motifs) {
        //basepair step
        if(m->mtype() == util::MotifType::HELIX) {
            motif = resources::Manager::instance().bp_step(m->end_ids()[0]);
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
    
    return 0;
}