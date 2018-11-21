//
// Created by Joseph Yesselman on 12/18/17.
//

#include "general_helix_sampler/general_helix_sampler.h"
#include "base/backtrace.hpp"
#include "util/cartesian_product.h"
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
    add_option("pdb", "", OptionType::STRING, false);
    add_option("motif", "", OptionType::STRING, false);
    add_option("start_bp", "", OptionType::STRING, true);
    add_option("end_bp", "", OptionType::STRING, true);
    add_option("seq", "", OptionType::STRING, true);
    add_option("all", false, OptionType::BOOL, false);
    add_option("get_ideal", false, OptionType::BOOL, false);
    add_option("get_dist_from_ideal", false, OptionType::BOOL, false);
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

    if(get_string_option("pdb") == "" && get_string_option("motif") == "") {
        throw GeneralHelixSamplerException("must supply pdb or motif");
    }

    auto start_bp = BasepairOP(nullptr);
    auto end_bp = BasepairOP(nullptr);
    auto start_bp_str = get_string_option("start_bp");
    auto end_bp_str = get_string_option("end_bp");

    if(get_string_option("pdb") != "") {
        auto pdb = get_string_option("pdb");
        auto rs = RM::instance().get_structure(pdb, "pdb");
        start_bp = rs->get_basepair(start_bp_str)[0];
        end_bp = rs->get_basepair(end_bp_str)[0];
    }
    else {
        auto motif_file = get_string_option("motif");
        auto m = file_to_motif(motif_file);
        start_bp = m->get_basepair(start_bp_str)[0];
        end_bp = m->get_basepair(end_bp_str)[0];
    }


    start_bp->bp_type("cW-W");
    end_bp->bp_type("cW-W");

    auto mf = MotifFactory();
    auto start = mf.motif_from_bps(BasepairOPs{start_bp, end_bp});
    start->ends(BasepairOPs{start_bp, end_bp});
    auto ref_m = mf.ref_motif();

    start->name("start");
    start->block_end_add(-1);

    RM::instance().register_motif(start);

    auto seq = get_string_option("seq");
    auto struc = _generate_structure(seq);

    auto bp_steps = get_motifs_from_seq_and_ss(seq, struc);

    auto mt = std::make_shared<MotifTree>();
    mt->set_option_value("sterics", false);

    mt->add_motif(start);
    mt->add_motif(bp_steps[0]);
    auto dist1 = start->ends()[0]->d().distance(start->ends()[1]->d());
    auto dist2 = mt->get_node(1)->data()->ends()[1]->d().distance(start->ends()[1]->d());

    if(dist1 < dist2) {
        mt->get_node(0)->data()->ends()[0]->flip();
    }

    auto new_start = mt->get_node(0)->data();

    if(! get_bool_option("all")) {
        std::cout << _get_hit_count(new_start, seq, struc) << std::endl;
        return;
    }

    auto pairs = std::vector<Strings>();
    pairs.push_back(Strings{"A", "U"});
    pairs.push_back(Strings{"U", "A"});
    pairs.push_back(Strings{"C", "G"});
    pairs.push_back(Strings{"G", "C"});

    auto all_pairs = std::vector<std::vector<Strings>>((int)(seq.size() / 2));
    for(int i = 0; i < all_pairs.size(); i++) {
        all_pairs[i] = pairs;
    }

    auto pair_iterator = CartesianProduct<Strings>(all_pairs);
    auto current = std::vector<Strings>();
    auto new_seq = String();
    auto seq1 = String();
    auto seq2 = String();

    while (!pair_iterator.end()) {
        current = pair_iterator.next();
        seq1 = "";
        seq2 = "";
        for (auto const & p : current) {
            seq1 += p[0];
            seq2 = p[1] + seq2;
        }
        new_seq = seq1 + "&" + seq2;
        std::cout << new_seq << " " <<  _get_hit_count(new_start, new_seq, struc) << std::endl;
    }

}

int
GeneralHelixSampler::_get_hit_count(
        MotifOP const & start,
        String const & seq,
        String const & ss) {

    auto bp_steps = get_motifs_from_seq_and_ss(seq, ss);
    auto mt = std::make_shared<MotifTree>();
    mt->set_option_value("sterics", false);
    mt->add_motif(start);
    int i = 0;
    for(auto const & e : start->ends()) {
        std::cout << e->name() << std::endl;
    }
    for(auto const & m : bp_steps) {
        mt->add_motif(m);
    }

    if(get_bool_option("get_ideal")) {
        auto mt2 = std::make_shared<MotifTree>();
        for(auto const & n : *mt) {
            if(n->data()->name() == "start") { continue; }
            mt2->add_motif(n->data());
        }
        mt2->to_pdb("ideal.pdb", 1, 1);
        exit(0);
    }

    if(get_bool_option("get_dist_from_ideal")) {
        auto target = mt->get_node(0)->data()->ends()[1];
        auto current = mt->last_node()->data()->ends()[1];

        std::cout << target->d().to_str() << "," << current->d().to_str() << ",";
        std::cout << target->d().distance(current->d()) << ",";
        std::cout << target->r().to_str() << "," << current->r().to_str() << ",";
        std::cout << target->r().difference(current->r()) << ",";
        std::cout << target->d().distance(current->d()) + 2*target->r().difference(current->r()) << std::endl;
        exit(0);

    }

    auto mset = std::make_shared<MotifStateEnsembleTree>(mt);
    tfs_.setup(mset, 0, mt->last_node()->index(), 1, 1);
    return tfs_.run();

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