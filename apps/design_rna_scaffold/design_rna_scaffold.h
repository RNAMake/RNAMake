//
// Created by Joseph Yesselman on 3/9/19.
//

#ifndef TEST_DESIGN_RNA_SCAFFOLD_H
#define TEST_DESIGN_RNA_SCAFFOLD_H

#include <stdio.h>
#include <regex>

#include <data_structure/graph_base.h>
#include "base/application.hpp"
#include "util/steric_lookup.hpp"
#include "util/visit_struct.h"
#include "util/visit_struct_intrusive.h"
#include "motif_search/search.h"
#include "motif_search/solution_topology.h"
#include <motif_search/solution_filter.h>
#include "motif_data_structure/motif_graph.h"
#include "sequence_optimization/sequence_optimizer_3d.hpp"
#include <thermo_fluctuation/graph/simulation.h>

#include <CLI/CLI.hpp>

struct GraphIndexes {
    data_structure::NodeIndexandEdge start = data_structure::NodeIndexandEdge();
    data_structure::NodeIndexandEdge end = data_structure::NodeIndexandEdge();
};

class DesignRNAScaffold : public base::Application {
public:

    DesignRNAScaffold();

public: // application functions

    void
    setup_options() override ;

    void
    parse_command_line(
            int,
            const char **) override ;

    void
    run() override ;

public: // getters

//TODO: Figure out a way to include this back

//    base::LogLevel
//    log_level() const { //added by CJ
//        return base::log_level_from_str(parameters_.core.log_level);
//    }

private: // setup functions

    void
    setup();

    void
    _setup_from_pdb();

    motif_search::SearchOP
    _setup_search();

    motif_search::ProblemOP
    _setup_problem();

    motif_search::SolutionTopologyTemplateOP
    _setup_sol_template_from_path(
            String const &);

    motif_search::SolutionFilterOP
    _setup_sol_filter(
            String const &);

    void
    _check_bp(
            String const &,
            structure::RNAStructureOP const &,
            String const &) const;

    void
    _setup_extra_pdbs();


private: // run functions

    motif_data_structure::MotifGraphOP
    _get_motif_graph_solution();

    void
    _get_graph_indexes_after_bp_steps(
            motif_data_structure::MotifGraph const &,
            GraphIndexes const &,
            GraphIndexes & /* return */);

    motif_data_structure::MotifGraphOP
    _perform_sequence_opt(
            motif_data_structure::MotifGraph const &,
            GraphIndexes const &);

    motif_data_structure::MotifGraphOP
    _perform_thermo_fluc_sim(
            motif_data_structure::MotifGraph &,
            GraphIndexes const &);

    void
    _record_solution(
            motif_data_structure::MotifGraph &);

    void
    _fix_flex_helices_mtype(
            motif_data_structure::MotifGraphOP);

    void
    _get_motif_names(
            motif_data_structure::MotifGraphOP);

    void
    _build_new_ensembles(
            String const &);

    void
    _get_mseg(
            motif_data_structure::MotifGraph &,
            motif_data_structure::MotifStateEnsembleGraph & /* return */,
            std::map<int, int> & /* return */);

private:
    struct SolutionInfo {
        // from motif search
        int design_num = 0;
        float design_score = 0.0f;
        String designable_sequence = "";
        String dot_bracket = "";
        String motif_names = "";
        // from sequence opt
        int seqeunce_opt_num = 0;
        String sequence = "";
        float sequence_opt_score = 0.0f;
        // from thermo fluc opt
        float thermo_fluc_best_score = 0.0f;
        int thermo_fluc_hits = 0.0;

        // listed in order of appearance
        static
        String
        col_names() {
            auto cols = String("");
            cols += "design_num,design_score,design_sequence,design_structure,motifs_uses,opt_num,";
            cols += "opt_sequence,opt_score,thermo_fluc_best_score,hit_count";
            return cols;
        }

        friend
        std::ostream&
        operator<<( std::ostream &output, SolutionInfo const & si) {
            output << si.design_num << "," << si.design_score << "," << si.designable_sequence << ",";
            output << si.dot_bracket << "," << si.motif_names << "," << si.seqeunce_opt_num << ",";
            output << si.sequence << "," << si.sequence_opt_score << "," << si.thermo_fluc_best_score;
            output << "," << si.thermo_fluc_hits;
            return output;
        }
    };


public:
    CLI::App app_; // added by CJ 08/20. has to be public to get the --help to work

private:

    // must be initialized a runtime
    resources::Manager & rm_;
    // general vars
    GraphIndexes starting_indexes_ = GraphIndexes();
    motif_data_structure::MotifStateGraphOP msg_;
    motif_data_structure::MotifGraphOP mg_;
    SolutionInfo sol_info_ = SolutionInfo();
    util::StericLookupNewOP lookup_;
    // search vars
    motif_search::SearchOP search_;
    motif_search::ProblemOP problem_;

    // sequence opt vars
    sequence_optimization::SequenceOptimizer3DOP seq_optimizer_;

    // thermo sim vars
    thermo_fluctuation::graph::SimulationOP thermo_sim_;

    // other vars
    std::ofstream out_, score_out_;
    std::map<String, motif::MotifStateEnsembleOP> new_motif_ensembles_;
private:

  struct Parameters;
};

#endif //TEST_DESIGN_RNA_SCAFFOLD_H