//
// Created by Joseph Yesselman on 4/12/18.
//

#include "base/backtrace.hpp"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_state_graph.hpp"
#include "motif_state_search/motif_state_monte_carlo.h"
#include "sequence_optimizer/sequence_optimizer_3d.hpp"
#include "apt_stablization/apt_stablization.h"



APTStablization::APTStablization() : Application() {}

void
APTStablization::setup_options() {
    add_option("pdb", String(""), OptionType::STRING, true);

    // general options
    add_option("out_file", "default.out", OptionType::STRING, false);
    add_option("score_file", "default.scores", OptionType::STRING, false);
    add_option("designs", 1, OptionType::INT, false);
    add_option("only_existing_motifs", false, OptionType::BOOL, false);
    add_option("only_ideal", false, OptionType::BOOL, false);
}

void
APTStablization::parse_command_line(
        int argc,
        const char **argv) {

    Application::parse_command_line(argc, argv);
}


void
APTStablization::run() {

    // add motif to resource manager
    RM::instance().add_motif(get_string_option("pdb"), "aptamer", MotifType::TWOWAY);
    std::cout << "APT STABLIZATION: loaded pdb from file: " << get_string_option("pdb") << std::endl;

    auto m = RM::instance().motif("aptamer");
    if(m->ends().size() < 2) {
        throw APTStablizationException("aptamer does not have not have two basepair ends");
    }

    //auto start_path = String("flex_helices,aptamer,flex_helices,twoway,flex_helices,twoway,flex_helices");
    auto start_path = String("flex_helices,twoway,flex_helices,aptamer,flex_helices,twoway,flex_helices");
    if(get_bool_option("only_existing_motifs")) {
        start_path = "flex_helices,aptamer,flex_helices,existing,flex_helices,existing,flex_helices";
    }

    auto ms_libraries = _get_libraries(start_path);

    // setup initial graph
    auto ttr = RM::instance().motif("GAAA_tetraloop", "", "A229-A245");
    auto bp_step_1 = RM::instance().bp_step("GC_LL_GC_RR");
    auto bp_step_2 = RM::instance().bp_step("CG_LL_CG_RR");
    auto msg = std::make_shared<MotifStateGraph>();
    msg->add_state(bp_step_1->get_state());
    msg->add_state(bp_step_2->get_state());
    msg->add_state(ttr->get_state());

    _setup_sterics(msg);

    auto start_end_pos = msg->get_node(2)->data()->get_end_index("A222-A251");
    auto end_end_pos = msg->get_node(2)->data()->get_end_index("A149-A154");

    bool target_an_aligned_end = false;
    if(end_end_pos == msg->get_node(2)->data()->block_end_add()) {
        target_an_aligned_end = true;
    }

    std::ofstream out, sf_out;
    sf_out.open(get_string_option("score_file"));
    sf_out << "design_num,design_score,design_sequence,design_structure,motifs_uses,opt_num,";
    sf_out << "opt_sequence,opt_score,eterna_score" << std::endl;

    out.open(get_string_option("out_file"));

    auto mc = MotifStateMonteCarlo(ms_libraries);
    mc.setup(msg, 2, 2, start_end_pos, end_end_pos, target_an_aligned_end);
    mc.lookup(lookup_);
    mc.start();
    mc.set_option_value("max_solutions", 10000);


    auto designs = 0;

    while(! mc.finished()) {
        auto sol = mc.next();
        if(sol == nullptr) {
            std::cout << "no solutions" << std::endl;
            exit(0);
        }

        auto mg = sol->mg;
        auto motif_names = String("");
        for(auto & n : *mg) {
            if(n->data()->name().substr(0,5) == "HELIX") {
                n->data()->mtype(MotifType::HELIX);
            }
            if(n->data()->name().substr(0,4) == "GAAA") {
                n->data()->mtype(MotifType::NWAY);
            }

            motif_names += n->data()->name() + ";";
        }
        mg->replace_ideal_helices();

        auto end_node = mg->get_node(2);
        auto target_state = end_node->data()->ends()[end_end_pos]->state();
        auto end_i = end_end_pos;
        auto partner = end_node->connections()[end_i]->partner(end_node->index());
        auto scorer = std::make_shared<InternalTargetScorer>(2, end_end_pos, partner->index(), 1, target_an_aligned_end);

        if(get_bool_option("only_ideal")) {
            sf_out << designs << "," << sol->score << "," << mg->designable_sequence() << "," ;
            sf_out << mg->dot_bracket() << "," << motif_names << ",";
            sf_out << std::endl;
            designs += 1;
            out << mg->to_str() << std::endl;
            mg->to_pdb("design." + std::to_string(designs) + ".pdb", 1, 1, 1);

            if(designs >= get_int_option("designs")) {
                std::cout << "Found " << designs << " designs, Finished!" << std::endl;
                break;
            }

            continue;
        }

        auto optimizer = SequenceOptimizer3D();
        optimizer.set_option_value("verbose", true);
        optimizer.set_option_value("cutoff", 7.0f);

        auto sols = optimizer.get_optimized_sequences(mg, scorer);
        if(sols[0]->dist_score > 7) {
            continue;
        }

        auto copy_mg = std::make_shared<MotifGraph>(*mg);
        auto dss = copy_mg->designable_secondary_structure();
        dss->replace_sequence(sols[0]->sequence);
        copy_mg->replace_helical_sequence(dss);
        copy_mg->to_pdb("design." + std::to_string(designs) + ".pdb", 1, 1, 1);

        sf_out << designs << "," << sol->score << "," << copy_mg->designable_sequence() << "," ;
        sf_out << copy_mg->dot_bracket() << "," << motif_names << ",";
        sf_out << 1 << "," << sols[0]->sequence << "," << sols[0]->dist_score << "," << sols[0]->eterna_score;
        sf_out << std::endl;
        out << copy_mg->to_str() << std::endl;
        designs += 1;

        if(designs >= get_int_option("designs")) {
            std::cout << "Found " << designs << " designs, Finished!" << std::endl;
            break;
        }


    }

    sf_out.close();
    out.close();

}


std::vector<MotifStateOPs>
APTStablization::_get_libraries(
        String const & motif_path) {
    auto spl = base::split_str_by_delimiter(motif_path, ",");
    auto i = 0;
    auto libraries = std::vector<MotifStateOPs>();
    auto motif_states = MotifStateOPs();
    for(auto const & name : spl) {
        if(name.length() < 2) { continue; }
        if(name == "ideal_helices_min" || name == "unique_twoway" || name == "tcontact" ||
           name == "twoway" || name == "flex_helices" || name == "existing") {
            auto ms_lib =  MotifStateSqliteLibrary(name);
            ms_lib.load_all();
            motif_states = MotifStateOPs();
            if(name == "flex_helices") {
                for (auto const & ms : ms_lib) {
                    if(ms->size() < 20) {
                        motif_states.push_back(ms);
                    }
                }

            }
            else {
                for (auto const & ms : ms_lib) { motif_states.push_back(ms); }
            }
            libraries.push_back(motif_states);
        }
        else {
            auto m = RM::instance().motif(name);
            motif_states = MotifStateOPs();

            for(auto const & end : m->ends()) {
                try {
                    auto m_new = RM::instance().motif(name, "", end->name());
                    m_new->new_res_uuids();
                    motif_states.push_back(m_new->get_state());
                }
                catch (...) {}
            }
            if(motif_states.size() == 0) {
                throw APTStablizationException("no viable aptamer conformations");
            }

            libraries.push_back(motif_states);
        }
    }
    return libraries;
}

void
APTStablization::_setup_sterics(
        MotifStateGraphOP msg) {
    auto beads = Points();
    for (auto & n : *msg) {
        for(auto const & b : n->data()->cur_state->beads()) {
            beads.push_back(b);
        }
    }
    lookup_.add_points(beads);
}

int main(int argc, const char *argv[]) {
    //must add this for all apps!
    std::set_terminate(base::print_backtrace);

    String base_path = base::base_dir() + "/rnamake/lib/RNAMake/apps/simulate_tectos/resources/";
    RM::instance().add_motif(base_path+"GAAA_tetraloop");

    auto app = APTStablization();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;

}
