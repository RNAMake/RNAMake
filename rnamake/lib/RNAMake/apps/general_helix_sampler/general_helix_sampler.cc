//
// Created by Joseph Yesselman on 8/30/17.
//

#include "base/backtrace.hpp"
#include "secondary_structure/secondary_structure_parser.h"
#include "motif_data_structures/motif_tree.h"
#include "resources/resource_manager.h"
#include "general_helix_sampler/general_helix_sampler.h"

GeneralHelixSampler::GeneralHelixSampler() : Application(),
        tfs_(ThermoFlucSimulation())
{}

void
GeneralHelixSampler::setup_options() {
    add_option("seq", "", OptionType::STRING);
    add_option("pdb", "", OptionType::STRING);
    add_option("start_bp", "", OptionType::STRING);
    add_option("end_bp", "", OptionType::STRING);

    add_cl_options(tfs_.options(), "simulation");
}

void
GeneralHelixSampler::parse_command_line(
        int argc,
        const char ** argv) {

    Application::parse_command_line(argc, argv);
    cl_parser_.assign_options(cl_options_, tfs_.options(), "simulation");
    tfs_.update_var_options();
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
            motif = RM::instance().motif("", m->end_ids()[0]);
            motifs.push_back(motif);
        }
        else {
            throw std::runtime_error("not supported");
        }
    }

    return motifs;
}



void
GeneralHelixSampler::run() {

    auto pdb = get_string_option("pdb");
    auto m = RM::instance().get_structure(get_string_option("pdb"), "scaffold");

    auto start_bp = m->get_basepair(get_string_option("start_bp"))[0];
    start_bp->flip();
    auto end_bp   = m->get_basepair(get_string_option("end_bp"))[0];
    end_bp->bp_type("cW-W");

    auto seq = get_string_option("seq");
    auto length = (seq.length()-1)/2;
    auto ss = String("");
    for(int i = 0; i < length; i++) { ss += "("; }
    ss += "&";
    for(int i = 0; i < length; i++) { ss += ")"; }

    auto mf = MotifFactory();
    auto start_motif = mf.motif_from_bps(BasepairOPs{start_bp, end_bp});
    start_motif->block_end_add(-1);

    RM::instance().register_motif(start_motif);

    auto mt = std::make_shared<MotifTree>();
    mt->set_option_value("sterics", false);
    mt->add_motif(start_motif);
    auto bp_steps = get_motifs_from_seq_and_ss(seq, ss);
    for(auto const & bp_step : bp_steps) {
        mt->add_motif(bp_step);
    }
    auto mset = std::make_shared<MotifStateEnsembleTree>(mt);
    tfs_.setup(mset, 0, mset->last_node()->index(), 1, 1);
    //tfs_.set_option_value("steps", 100000);
    auto count = tfs_.run();
    std::cout << count << std::endl;
    exit(0);



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


