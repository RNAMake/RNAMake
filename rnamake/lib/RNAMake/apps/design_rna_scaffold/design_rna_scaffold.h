//
// Created by Joseph Yesselman on 3/9/19.
//

#ifndef TEST_DESIGN_RNA_SCAFFOLD_H
#define TEST_DESIGN_RNA_SCAFFOLD_H

#include <stdio.h>
#include <data_structure/graph_base.h>
#include "base/application.hpp"
#include "util/steric_lookup.hpp"
#include "motif_search/search.h"
#include "motif_search/solution_topology.h"
#include <motif_search/solution_filter.h>
#include "motif_data_structure/motif_graph.h"
#include "sequence_optimization/sequence_optimizer_3d.hpp"


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
    _setup_sterics();

    void
    _setup_from_pdb();

    void
    _setup_from_mg();

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
            structure::RNAStructureOP,
            String const &);

private:
    struct Parameters {
        String pdb, start_bp, end_bp, mg;
        String starting_helix, ending_helix, search_type, motif_path;
        String out_file, score_file, solution_filter, new_ensembles;
        bool skip_sequence_optimization, no_basepair_checks, no_mg_file;
        bool all_designs, dump_pdbs, dump_scaffold_pdbs;
        float search_cutoff;
        int search_max_size, designs;
        int max_helix_length, min_helix_length;
        float scaled_score_d, scaled_score_r;
        //scoring related parameters
        String exhaustive_scorer, mc_scorer;

    };

private:
    std::ofstream out_, score_out_;
    motif_search::SearchOP search_;
    motif_search::ProblemOP problem_;
    motif_data_structure::MotifStateGraphOP msg_;
    data_structure::NodeIndexandEdge start_, end_;
    sequence_optimization::SequenceOptimizer3D optimizer_;
    resources::Manager & rm_;
    Parameters parameters_;
    String motif_names_;
    std::map<String, motif::MotifStateEnsembleOP> new_motif_ensembles_;


};



#endif //TEST_DESIGN_RNA_SCAFFOLD_H
