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

SequenceOptimizationBenchmarks::SequenceOptimizationBenchmarks():
        parameters_(SequenceOptimizationParameters()){}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// app functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
SequenceOptimizationBenchmarks::setup_options() {
    add_option("problem", String("TTR"), OptionType::STRING, false);
    add_option("helices", String("ideal_helices"), OptionType::STRING, false);
    add_option("rounds", 1, OptionType::INT, false);
    add_option("motifs", 3, OptionType::INT, false);
    add_option("min_helix_size", 6, OptionType::INT, false);
    add_option("max_helix_size", 18, OptionType::INT, false);

    add_option("out_file", "results.csv", OptionType::STRING, false);

}

void
SequenceOptimizationBenchmarks::parse_command_line(
        int argc,
        const char **argv) {

    Application::parse_command_line(argc, argv);

    parameters_.problem        = get_string_option("problem");
    parameters_.helices        = get_string_option("helices");
    parameters_.rounds         = get_int_option("rounds");
    parameters_.motifs         = get_int_option("motifs");
    parameters_.min_helix_size = get_int_option("min_helix_size");
    parameters_.max_helix_size = get_int_option("max_helix_size");

}

void
SequenceOptimizationBenchmarks::run() {
    //auto problem_factory = std::make_shared<TTRProblemFactory>();
    auto problem_factory = std::make_shared<RibosomeTetherProblemFactory>();
    auto ms_libraries = _get_libraries();
    auto timer = Timer();

    std::ofstream f_out;
    f_out.open(get_string_option("out_file"));

    for(int i = 0; i < parameters_.rounds; i++) {
        auto problem = problem_factory->get_problem();
        auto search  = _get_search(problem, ms_libraries);

        timer.start();
        auto sol = search->next();
        auto time_1 = timer.end();

        if (sol == nullptr) { continue; }
        auto & motif_names = _get_motif_names(sol->mg);

        // fix flex helices
        for(auto & n : *sol->mg) {
            if(n->data()->name().substr(0,1) == "H") {
                n->data()->mtype(MotifType::HELIX);
            }
        }

        sol->mg->write_pdbs();
        exit(0);

        sol->mg->replace_ideal_helices();
        auto optimizer = _get_optimizer(problem, sol->mg);

        timer.start();
        auto sols = optimizer->get_optimized_sequences(sol->mg);
        auto time_2 = timer.end();

        f_out << sol->score << "," << time_1 << "," << sols[0]->dist_score << "," << time_2 << "," << motif_names << std::endl;
        f_out.flush();
    }

    f_out.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::vector<MotifStateOPs>
SequenceOptimizationBenchmarks::_get_libraries() {

    auto num = parameters_.motifs*2 + 1;
    auto motif_lib_names = Strings(num);
    auto motif_types = std::vector<MotifType>(num);
    for(int i = 0; i < num; i++) {
        if (i % 2 == 0) {
            motif_lib_names[i] = parameters_.helices;
            motif_types[i]     = MotifType::HELIX;
        }
        else            {
            motif_lib_names[i] = "twoway";
            motif_types[i]     = MotifType::TWOWAY;
        }
    }

    auto libraries = std::vector<MotifStateOPs>();
    auto motif_states = MotifStateOPs();
    int i = 0;
    for(auto const & name : motif_lib_names) {
        auto ms_lib =  MotifStateSqliteLibrary(name);
        ms_lib.load_all();
        motif_states = MotifStateOPs();

        if(motif_types[i] == MotifType::HELIX) {
            for (auto const & ms : ms_lib) {
                if(ms->size() >= parameters_.min_helix_size && ms->size() <= parameters_.max_helix_size) {

                    motif_states.push_back(ms);
                }
            }
        }
        else {
            for (auto const & ms : ms_lib) { motif_states.push_back(ms); }
        }
        libraries.push_back(motif_states);
        i++;
    }

    return libraries;
}


MotifStateMonteCarloOP
SequenceOptimizationBenchmarks::_get_search(
        SequenceOptProblemOP problem,
        std::vector<MotifStateOPs> const & ms_libraries) {

    auto mc = std::make_shared<MotifStateMonteCarlo>(ms_libraries);
    mc->setup(problem->msg, problem->start.ni, problem->end.ni,
              problem->start.ei, problem->end.ei, problem->target_an_aligned_end);
    mc->set_option_value("accept_score", 5.0f);
    mc->lookup(*problem->lookup);
    mc->start();

    return mc;
}

SequenceOptimizer3DOP
SequenceOptimizationBenchmarks::_get_optimizer(
        SequenceOptProblemOP problem,
        MotifGraphOP mg) {

    auto end_node = mg->get_node(problem->end.ni);
    auto target_state = end_node->data()->ends()[problem->end.ei]->state();
    auto partner = end_node->connections()[problem->end.ei]->partner(end_node->index());
    auto scorer = std::make_shared<InternalTargetScorer>(problem->end.ni, problem->end.ei, partner->index(), 1,
                                                         problem->target_an_aligned_end);

    auto optimizer = std::make_shared<SequenceOptimizer3D>();
    //optimizer.set_option_value("verbose", true);
    optimizer->set_option_value("cutoff", 0.01f);
    optimizer->set_option_value("return_lowest", true);
    optimizer->set_scorer(scorer);

    return optimizer;

}

String const &
SequenceOptimizationBenchmarks::_get_motif_names(
        MotifGraphOP mg) {
    motif_names_ = "";
    for (auto const & n : *mg) {
        motif_names_ += n->data()->name() + ";";
    }
    return motif_names_;
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