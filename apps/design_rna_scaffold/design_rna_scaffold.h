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
#include "motif_search/search.h"
#include "motif_search/solution_topology.h"
#include <motif_search/solution_filter.h>
#include "motif_data_structure/motif_graph.h"
#include "sequence_optimization/sequence_optimizer_3d.hpp"
#include <thermo_fluctuation/graph/simulation.h>

#include <CLI/CLI.hpp>

class DesignRNAScaffold : public base::Application {
public:

    DesignRNAScaffold();

public: // application functions

    void
    setup_options();

    void
    parse_command_line(
            int,
            const char **);

    void
    run();

private:

    void
    setup();

private:

    bool
    _get_motif_graph_solution();

    bool
    _get_sequence_optimization_solution(

        );

    void
    _setup_sterics();

    void
    _setup_from_pdb();

    std::vector<motif::MotifStateOPs>
    _get_libraries();

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
    _record_solution(
            motif_data_structure::MotifGraphOP,
            motif_search::SolutionOP,
            sequence_optimization::OptimizedSequenceOP,
            int,
            int,
            int);

    void
    _fix_flex_helices_mtype(
            motif_data_structure::MotifGraphOP);

    String const &
    _get_motif_names(
            motif_data_structure::MotifGraphOP);

    void
    _build_new_ensembles(
            String const &);

private:

    struct EnsembleConversionResults {
        inline
        EnsembleConversionResults(
                motif_data_structure::MotifStateEnsembleGraphOP n_mseg,
                std::map<int, int> const & n_index_hash):
                mseg(n_mseg),
                index_hash(n_index_hash) {}

        motif_data_structure::MotifStateEnsembleGraphOP mseg;
        std::map<int, int> index_hash;
    };

    typedef std::shared_ptr<EnsembleConversionResults> EnsembleConversionResultsOP;

    EnsembleConversionResultsOP
    _get_mseg(
            motif_data_structure::MotifGraphOP);


private:
    void
    check_bp(
            String const &,
            structure::RNAStructureOP const &,
            String const &);

private:


    struct Parameters {
        // required options will exit if not supplied
        struct Core {
            String pdb = "";
            String start_bp = "";
            String end_bp = "";
            int designs = 1;
        };
        // options related to what will be outputted
        struct IO {
            String out_file = "default.out";
            String score_file = "default.scores";
            String new_ensembles_file = "";
            bool dump_pdbs = false;
            bool dump_scaffold_pdbs = false;
            bool dump_intermediate_pdbs = false;
            bool no_out_file = false;
        // options related to how the initial motif search is performed
        };
        struct Search {
            String type = "path_finding";
            String motif_path = "";
            String starting_helix = "";
            String ending_helix = "";
            String solution_filter = "NotFilter";
            float cutoff = 5.0f;
            int max_helix_length = 99;
            int min_helix_length = 4;
            int max_size = 9999;
            bool no_basepair_checks = false;
            String exhaustive_scorer, mc_scorer;
            float scaled_score_d, scaled_score_r;

        };

        struct SequenceOpt {
            bool skip;
            bool sequences_per_design;

        };

        struct ThermoFluc {
            bool skip;

        };

        Core core = Core();
        IO io = IO();
        Search search = Search();
        SequenceOpt seq_opt = SequenceOpt();
        ThermoFluc thermo_fluc = ThermoFluc();
    };

public:
    CLI::App app_; // added by CJ 08/20. has to be public to get the --help to work

private:
    // must be initialized a runtime
    resources::Manager & rm_;
    // general vars
    Parameters parameters_ = Parameters();
    data_structure::NodeIndexandEdge start_ = data_structure::NodeIndexandEdge();
    data_structure::NodeIndexandEdge end_ = data_structure::NodeIndexandEdge();
    motif_data_structure::MotifStateGraphOP msg_;
    motif_data_structure::MotifGraphOP mg_;
    motif_data_structure::MotifGraphOP mg_w_sol_;


    // search vars
    motif_search::SearchOP search_;
    motif_search::ProblemOP problem_;
    motif_search::SolutionOP solution_;
    int design_num_ = 0;
    // sequence opt vars
    sequence_optimization::SequenceOptimizer3DOP seq_optimizer_;

    // thermo sim vars
    thermo_fluctuation::graph::SimulationOP thermo_sim_;
    // other vars
    std::ofstream out_, score_out_;
    String motif_names_ = "";
    std::map<String, motif::MotifStateEnsembleOP> new_motif_ensembles_;


};



#endif //TEST_DESIGN_RNA_SCAFFOLD_H
