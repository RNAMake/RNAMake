//
// Created by Joseph Yesselman on 12/18/17.
//

#include "general_helix_sampler/general_helix_sampler.h"
#include "base/backtrace.hpp"
#include "secondary_structure/secondary_structure_parser.h"
#include "motif_data_structures/motif_tree.h"
#include "motif_data_structures/motif_state_ensemble_tree.h"
#include "resources/resource_manager.h"
#include "thermo_fluctuation/thermo_fluc_sampler.h"

GeneralHelixSampler::GeneralHelixSampler() : Application(),
        tfs_(ThermoFlucSimulation())
{}

// application setups functions ////////////////////////////////////////////////////////////////////

void
GeneralHelixSampler::setup_options() {
    add_option("pdb", "", OptionType::STRING, true);
    add_option("start_bp", "", OptionType::STRING, true);
    add_option("end_bp", "", OptionType::STRING, true);
    add_option("seq", "", OptionType::STRING, true);

}

void
GeneralHelixSampler::parse_command_line(
        int argc,
        const char ** argv) {
    Application::parse_command_line(argc, argv);
}

String
GeneralHelixSampler::_generate_structure(
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
GeneralHelixSampler::get_motifs_from_seq_and_ss(
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



void
GeneralHelixSampler::run() {
    auto pdb = get_string_option("pdb");
    //auto structure = _generate_structure(sequence);
    auto rs =  RM::instance().get_structure(pdb, "pdb");
    auto start_bp_str = get_string_option("start_bp");
    auto end_bp_str = get_string_option("end_bp");

    auto start_bp = rs->get_basepair(start_bp_str)[0];
    auto end_bp = rs->get_basepair(end_bp_str)[0];

    start_bp->bp_type("cW-W");
    end_bp->bp_type("cW-W");

    auto mf = MotifFactory();
    auto start = mf.motif_from_bps(BasepairOPs{start_bp, end_bp});
    start->name("start");
    start->block_end_add(-1);

    RM::instance().register_motif(start);

    auto seq = get_string_option("seq");
    auto struc = _generate_structure(seq);

    auto bp_steps = get_motifs_from_seq_and_ss(seq, struc);

    auto mt = std::make_shared<MotifTree>();
    auto mt2 = std::make_shared<MotifTree>();
    mt->set_option_value("sterics", false);

    mt->add_motif(start);
    for(auto const & m : bp_steps) {
        mt->add_motif(m);
    }

    /*for(int i = 0; i < bp_steps.size(); i++) {
        mt2->add_motif(mt->get_node(i+1)->data());
    }

    mt2->to_pdb("merged.pdb", 1);
    exit(0);
    */
    auto mset = std::make_shared<MotifStateEnsembleTree>(mt);
    tfs_.setup(mset, 0, mt->last_node()->index(), 1, 1);
    std::cout << tfs_.run() << std::endl;

}




int main(int argc, const char * argv[]) {
    //must add this for all apps!
    std::set_terminate(print_backtrace);

    auto app = GeneralHelixSampler();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;


}