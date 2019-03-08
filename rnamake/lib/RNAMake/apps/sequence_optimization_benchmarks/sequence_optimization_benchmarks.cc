//
// Created by Joseph Yesselman on 3/7/19.
//

#include <chrono>

#include "base/backtrace.hpp"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_state_graph.hpp"
#include "motif_state_search/motif_state_monte_carlo.h"
#include "sequence_optimizer/sequence_optimizer_3d.hpp"
#include "sequence_optimization_benchmarks/sequence_optimization_benchmarks.h"

SequenceOptimizationBenchmarks::SequenceOptimizationBenchmarks() {}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// app functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
SequenceOptimizationBenchmarks::setup_options() {
    add_option("problem", String("TTR"), OptionType::STRING, false);
    add_option("helices", String("ideal"), OptionType::STRING, false);
    add_option("rounds", 1, OptionType::INT, false);

    add_option("out_file", "results.csv", OptionType::STRING, false);

}

void
SequenceOptimizationBenchmarks::parse_command_line(
        int argc,
        const char **argv) {

    Application::parse_command_line(argc, argv);
}

void
SequenceOptimizationBenchmarks::run() {
    auto msg = _get_starting_graph();
    auto lookup = StericLookup();
    _setup_sterics(msg, lookup);

    auto start_end_pos = msg->get_node(2)->data()->get_end_index("A222-A251");
    auto end_end_pos = msg->get_node(2)->data()->get_end_index("A149-A154");

    bool target_an_aligned_end = false;
    if(end_end_pos == msg->get_node(2)->data()->block_end_add()) {
        target_an_aligned_end = true;
    }

    for(int mode = 0; mode < 2; mode++) {

        auto motif_path = String();

        std::ofstream f_out;

        if(mode == 0) {
            f_out.open("results_ideal.dat");
            motif_path = "ideal_helices,twoway,ideal_helices,twoway,ideal_helices,twoway,ideal_helices";
        }
        else {
            f_out.open("results_avg.dat");
            motif_path = "avg_helices,twoway,avg_helices,twoway,avg_helices,twoway,avg_helices";
        }
        f_out << "design_score,design_time,opt_score,opt_time,motifs" << std::endl;
        auto ms_libraries = _get_libraries(motif_path);

        for (int i = 0; i < 100; i++) {
            msg = _get_starting_graph();
            auto mc = MotifStateMonteCarlo(ms_libraries);
            mc.setup(msg, 2, 2, start_end_pos, end_end_pos, target_an_aligned_end);
            mc.set_option_value("accept_score", 5.0f);
            mc.lookup(lookup);
            mc.start();

            auto start_1 = std::chrono::steady_clock::now();
            auto sol = mc.next();
            auto end_1 = std::chrono::steady_clock::now();
            auto time_1 = std::chrono::duration<double, std::milli> (end_1 - start_1).count();

            if(sol == nullptr) { continue; }

            auto mg = sol->mg;
            auto motif_names = String("");

            for(auto const & n : *mg) {
                motif_names += n->data()->name() + ";";
            }

            mg->replace_ideal_helices();

            auto end_node = mg->get_node(2);
            auto target_state = end_node->data()->ends()[end_end_pos]->state();
            auto end_i = end_end_pos;
            auto partner = end_node->connections()[end_i]->partner(end_node->index());
            auto scorer = std::make_shared<InternalTargetScorer>(2, end_end_pos, partner->index(), 1,
                                                                 target_an_aligned_end);

            auto optimizer = SequenceOptimizer3D();
            //optimizer.set_option_value("verbose", true);
            optimizer.set_option_value("cutoff", 0.01f);
            optimizer.set_option_value("return_lowest", true);

            auto start_2 = std::chrono::steady_clock::now();
            auto sols = optimizer.get_optimized_sequences(mg, scorer);
            auto end_2 = std::chrono::steady_clock::now();
            auto time_2 = std::chrono::duration<double, std::milli> (end_2 - start_2).count();
            f_out << sol->score << "," << time_1 << "," << sols[0]->dist_score << "," << time_2 << ",";
            f_out << motif_names << std::endl;
            f_out.flush();

        }
        f_out.close();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

MotifStateGraphOP
SequenceOptimizationBenchmarks::_get_starting_graph() {
    auto ttr = RM::instance().motif("GAAA_tetraloop", "", "A229-A245");
    auto bp_step_1 = RM::instance().bp_step("GC_LL_GC_RR");
    auto bp_step_2 = RM::instance().bp_step("CG_LL_CG_RR");
    auto msg = std::make_shared<MotifStateGraph>();
    msg->add_state(bp_step_1->get_state());
    msg->add_state(bp_step_2->get_state());
    msg->add_state(ttr->get_state());
    return msg;
}

void
SequenceOptimizationBenchmarks::_setup_sterics(
        MotifStateGraphOP msg,
        StericLookup & lookup) {
    auto beads = Points();
    for (auto & n : *msg) {
        for(auto const & b : n->data()->cur_state->beads()) {
            beads.push_back(b);
        }
    }
    lookup.add_points(beads);
}

std::vector<MotifStateOPs>
SequenceOptimizationBenchmarks::_get_libraries(
        String const & motif_path) {
    auto spl = split_str_by_delimiter(motif_path, ",");
    auto i = 0;
    auto libraries = std::vector<MotifStateOPs>();
    auto motif_states = MotifStateOPs();
    for(auto const & name : spl) {
        if(name.length() < 2) { continue; }
        if(name == "ideal_helices_min" || name == "unique_twoway" || name == "tcontact" || name == "ideal_helices" ||
           name == "twoway" || name == "flex_helices" || name == "existing" || name == "avg_helices") {
            auto ms_lib =  MotifStateSqliteLibrary(name);
            ms_lib.load_all();
            motif_states = MotifStateOPs();
            if(name == "ideal_helices") {
                for (auto const & ms : ms_lib) {
                    if(ms->size() > 4 && ms->size() < 18) {
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
                throw std::runtime_error("no viable aptamer conformations");
            }

            libraries.push_back(motif_states);
        }
    }
    return libraries;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int
main(
        int argc,
        const char ** argv) {

    //must add this for all apps!
    std::set_terminate(print_backtrace);

    //load extra motifs being used
    String base_path = base_dir() + "/rnamake/lib/RNAMake/apps/simulate_tectos/resources/";
    RM::instance().add_motif(base_path + "GAAA_tetraloop");

    auto app = SequenceOptimizationBenchmarks();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;

}