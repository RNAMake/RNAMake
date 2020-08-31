//
// Created by Joseph Yesselman on 3/9/19.
//

#include "base/backtrace.h"
#include "base/log.h"
#include "design_rna_scaffold/design_rna_scaffold.h"

#include <motif_search/solution_topology.h>
#include <motif_search/path_finding/search.h>
#include <motif_search/exhaustive/search.h>
#include <motif_search/monte_carlo/search.h>

#include <sequence_optimization/sequence_optimizer.h>
#include <thermo_fluctuation/graph/simulation.h>

DesignRNAScaffold::DesignRNAScaffold():
        base::Application(),
        rm_(resources::Manager::instance()),
        parameters_(Parameters()),
        app_("DesignRNAScaffold")
        {


}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// app functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
DesignRNAScaffold::setup_options() {

    app_.add_option_group("Core Inputs");
    app_.add_option_group("File Options");
    app_.add_option_group("Search Parameters");
    app_.add_option_group("Scoring Paramteres");
    app_.add_option_group("TBD");

    app_.add_option("--pdb",
                    parameters_.pdb,
                    "path to a PDB file with input RNA structure")
                    ->required()
                    ->check(CLI::ExistingFile)
                    ->group("Core Inputs");

    app_.add_option("--start_bp",
                        parameters_.start_bp,
                        "starting basepair to be used in structure format: [CHAIN ID][NT1 NUM] - [CHAIN ID][NT2 NUM]")
                    ->required()
                    ->group("Core Inputs");

    app_.add_option("--end_bp",
                    parameters_.end_bp,
                    "ending basepair to be used in structure format: [CHAIN ID][NT1 NUM] - [CHAIN ID][NT2 NUM]")
                    ->required()
                    ->group("Core Inputs");

    app_.add_option("--designs",
                        parameters_.designs,
                        "number of designs to create. Default is 1")
                    ->default_val(1)
                    ->check(CLI::PositiveNumber)
                    ->group("Core Inputs");

    app_.add_option("--out_file",
                    parameters_.out_file,
                    "output file for design(s)")
                    ->default_val("default.out")
                    ->group("File Options");

    parameters_.mg = "";

    app_.add_flag("--dump_pdbs",
                    parameters_.dump_pdbs,
                    "flag to dump intermediate pdbs TODO")
                    ->group("File Options");
                parameters_.dump_pdbs = false;

    app_.add_flag("--dump_scaffold_pdbs",
                    parameters_.dump_scaffold_pdbs,
                    "flag to output intermediate scaffold pdbs in the current directory")
                    ->group("File Options");
                    parameters_.dump_scaffold_pdbs = false;

    app_.add_option("--search_type",
                    parameters_.search_type,
                    "search type for traversing motif space")
                    ->default_val("path_finding")
                    ->check(CLI::IsMember(std::set<String>{"path_finding", "exhaustive", "mc"}))
                    ->group("Search Parameters");

    app_.add_option("--search_cutoff",
                    parameters_.search_cutoff,
                    "TODO")
                    ->default_val(5.0f)
                    ->group("Search Parameters");

    app_.add_option("--search_max_size",
                    parameters_.search_max_size,
                    "maximum number of steps for a design search")
                    ->default_val(999999)
                    ->check(CLI::PositiveNumber) //TODO put a limit to this?? CJ 08/20
                    ->group("Search Parameters");

    app_.add_flag("--skip_sequence_optimization",
                    parameters_.skip_sequence_optimization,
                    "flag to skip sequence optimization of the design")
                    ->group("Search Parameters");
                    parameters_.skip_sequence_optimization = false;

    app_.add_option("--score_file",
                    parameters_.score_file,
                    "name of output file containining scoring information for design")
                    ->default_val("default.scores")
                    ->group("File Options");


    app_.add_option("--solution_filter",
                    parameters_.solution_filter,
                    "TODO")
                    ->default_val("NotFilter")
                    ->check(CLI::IsMember(std::set<String>{"NoFilter","RemoveDuplicateHelices"}))
                    ->group("Search Parameters");

    // less common options
    app_.add_flag("--no_basepair_checks",
                    parameters_.no_basepair_checks,
                    "flag to disable basepair checks")
                    ->group("Search Parameters");
                    parameters_.no_basepair_checks = false;

    app_.add_option("--max_helix_length",
                parameters_.max_helix_length,
                    "maximum number of basepairs in a solution helix")
                    ->default_val(99)
                    ->group("Search Parameters");

    app_.add_option("--min_helix_length",
                    parameters_.min_helix_length,
                    "minimum number of basepairs in a solution helix")
                    ->default_val(4)
                    ->group("Search Parameters");

    app_.add_option("--starting_helix",
                    parameters_.starting_helix,
                    "starting helix for design solution. Format = [TODO]")
                    ->default_val("")
                    ->group("Search Parameters");

    app_.add_option("--ending_helix",
                    parameters_.ending_helix,
                    "ending helix for design solution. Format = [TODO]")
                    ->default_val("")
                    ->group("Search Parameters");

    app_.add_flag("--all_designs",
                    parameters_.all_designs,
                    "TBD")
                    ->default_val(false)
                    ->group("TBD");

    app_.add_option("--motif_path",
                    parameters_.motif_path,
                    "TBD")
                    ->default_val("")
                    ->group("TBD");

    app_.add_option("--new_ensembles",
                    parameters_.new_ensembles,
                    "TBD")
                    ->default_val("")
                    ->group("TBD");

    app_.add_flag("--no_mg_file",
                  parameters_.no_mg_file,
                  "TBD")
                  ->default_val(false)
                  ->group("TBD");

    app_.add_option("--exhaustive_scorer",
                    parameters_.exhaustive_scorer,
                    "TODO")
                    ->default_val("default")
                    ->group("Scoring Parameters");

    app_.add_option("--mc_scorer",
                    parameters_.mc_scorer,
                    "TODO")
                    ->default_val("default")
                    ->check(CLI::IsMember(std::set<String>{"default","scaled_scorer"}))
                    ->group("Scoring Parameters");

    app_.add_option("--scaled_score_d",
                    parameters_.scaled_score_d,
                    "TODO")
                    ->default_val(1.0f)
                    ->group("Scoring Parameters");

    app_.add_option("--scaled_score_r",
                    parameters_.scaled_score_r,
                    "TODO")
                    ->default_val(2.0f)
                    ->group("Scoring Parameters");
}

void
DesignRNAScaffold::parse_command_line(
        int argc,
        const char ** argv) {

        base::Application::parse_command_line(argc, argv);

    // core inputs
    // parameters_.pdb       = get_string_option("pdb");
    // parameters_.start_bp   = get_string_option("start_bp");
    // parameters_.end_bp    = get_string_option("end_bp");
    parameters_.mg         = get_string_option("mg");
    // common options
    parameters_.designs   = get_int_option("designs"); parameters_.dump_pdbs = get_bool_option("dump_pdbs"); parameters_.dump_scaffold_pdbs = get_bool_option("dump_scaffold_pdbs");
    parameters_.out_file = get_string_option("out_file"); parameters_.score_file = get_string_option("score_file");
    parameters_.solution_filter = get_string_option("solution_filter");
    parameters_.search_type = get_string_option("search_type");
    parameters_.motif_path = get_string_option("motif_path");
    parameters_.all_designs = get_bool_option("all_designs");
    parameters_.skip_sequence_optimization = get_bool_option("skip_sequence_optimization");
    // less common options
    parameters_.no_basepair_checks         = get_bool_option("no_basepair_checks");
    parameters_.no_mg_file = get_bool_option("no_mg_file");
    parameters_.starting_helix = get_string_option("starting_helix");
    parameters_.ending_helix = get_string_option("ending_helix");
    parameters_.search_cutoff = get_float_option("search_cutoff");
    parameters_.search_max_size = get_int_option("search_max_size");
    parameters_.new_ensembles = get_string_option("new_ensembles");
    parameters_.max_helix_length = get_int_option("max_helix_length");
    parameters_.min_helix_length = get_int_option("min_helix_length");

    // scoring related options
    parameters_.exhaustive_scorer = get_string_option("exhaustive_scorer");
    parameters_.mc_scorer = get_string_option("mc_scorer");
    parameters_.scaled_score_d = get_float_option("scaled_score_d");
    parameters_.scaled_score_r = get_float_option("scaled_score_r");

 }

void
DesignRNAScaffold::run() {
    if(parameters_.new_ensembles != "") { _build_new_ensembles(parameters_.new_ensembles); }

    if     (parameters_.pdb != "") { _setup_from_pdb(); }
    else if(parameters_.mg  != "") {}
    else                           {}

    out_.open(parameters_.out_file);
    score_out_.open(parameters_.score_file);

    score_out_ << "design_num,design_score,design_sequence,design_structure,motifs_uses,opt_num,";
    score_out_ << "opt_sequence,opt_score,eterna_score,hit_count" << std::endl;

    search_ = _setup_search();
    problem_ = _setup_problem();

    //sequence optimziation setup
    auto seq_optimizer = std::make_shared<sequence_optimization::SequenceOptimizer3D>();
    seq_optimizer->set_option_value("steps", 1);

    //thermo sim setup
    auto thermo_scorer = std::make_shared<thermo_fluctuation::graph::FrameScorer>();
    auto sterics = std::make_shared<thermo_fluctuation::graph::sterics::NoSterics>();
    auto sim = std::make_shared<thermo_fluctuation::graph::Simulation>(thermo_scorer, sterics);

    auto mg = msg_->to_motif_graph();
    auto sol_mg = motif_data_structure::MotifGraphOP(nullptr);
    mg->set_option_value("sterics", false);
    mg->increase_level();

    auto sol = motif_search::SolutionOP(nullptr);
    search_->setup(problem_);
    int i = 0;
    while(!search_->finished()) {

        try {
            sol = search_->next();

            if (sol == nullptr) { break; }

            sol_mg = sol->graph->to_motif_graph();
            if (sol_mg->size() == 0) { break; }
            _get_motif_names(sol_mg);

            if (i == 0) {
                LOG_INFO << "found a solution: " << motif_names_;
            } else if (i % 10 == 0) {
                LOG_INFO << "found " << i << " solutions ";
            }

            mg->add_motif_graph(*sol_mg, start_.node_index, start_.edge_index);
            auto mg_copy = std::make_shared<motif_data_structure::MotifGraph>(*mg);
            mg_copy->add_connection(mg->last_node()->index(), end_.node_index,
                                    mg->last_node()->data()->end_name(1),
                                    mg->get_node(end_.node_index)->data()->end_name(end_.edge_index));
            _fix_flex_helices_mtype(mg_copy);
            mg_copy->replace_ideal_helices();

            if(parameters_.skip_sequence_optimization) {
                _record_solution(mg_copy, sol, nullptr, 0, i, 0);
                mg->remove_level(1);
                continue;
            }

            auto last_m = mg_copy->get_node(end_.node_index)->connections()[end_.edge_index]->partner(end_.node_index);
            auto new_start = data_structure::NodeIndexandEdge{last_m->index(), 1};

            for (int j = 0; j < 1; j++) {

                // sequence optimization
                auto opt_seq_scorer = std::make_shared<sequence_optimization::InternalTargetScorer>(
                        new_start.node_index, new_start.edge_index, end_.node_index, end_.edge_index,
                        problem_->target_an_aligned_end);

                auto sols = seq_optimizer->get_optimized_sequences(mg_copy, opt_seq_scorer);
                mg_copy->replace_helical_sequence(sols[0]->sequence);

                // thermo sim
                auto r = _get_mseg(mg_copy);

                auto start = data_structure::NodeIndexandEdge{r->index_hash[new_start.node_index],
                                                              new_start.edge_index};
                auto end = data_structure::NodeIndexandEdge{r->index_hash[end_.node_index], end_.edge_index};

                sim->setup(*r->mseg, start, end);
                auto count = 0;
                for (int s = 0; s < 1000000; s++) {
                    count += sim->next();
                }
                _record_solution(mg_copy, sol, sols[0], count, i, j);

            }


            mg->remove_level(1);
        }
        catch(...) {
            mg->remove_level(1);
            continue;
        }

        i++;
        if(parameters_.designs <= i) {
            LOG_INFO << "found " << i << " designs, if you want more please use -designs num_of_design";
            exit(0);
        }
    }

    if(i == 0) {
        LOG_INFO << "no solutions found";
    }
    else {
        LOG_INFO << "no more solutions found";
        LOG_INFO << "found a total of " << i << " solution(s)";
    }

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
DesignRNAScaffold::_setup_from_pdb() {
    auto struc = rm_.get_structure(parameters_.pdb, "scaffold");
    LOGI << "loaded pdb from file: " << parameters_.pdb;

    if (parameters_.start_bp == "" || parameters_.end_bp == "") {
        LOG_ERROR << "must supply the name of the start_bp and end_bp using option -start_bp and "
                    "-end_bp respectively when using -pdb option";
    }

    check_bp(parameters_.start_bp, struc, "start");
    check_bp(parameters_.end_bp,  struc, "end");

    //auto m = std::make_shared<motif::Motif>()

    // TODO allow for construction of motif from RNA structure to allow for changes in base pair types
    rm_.add_motif(parameters_.pdb, "scaffold", util::MotifType::TCONTACT);
    auto m = rm_.motif("scaffold", "", parameters_.end_bp);

    auto ei1 = m->get_end_index(parameters_.start_bp);
    auto ei2 = m->get_end_index(parameters_.end_bp);

    start_ = data_structure::NodeIndexandEdge{0, ei1};
    end_   = data_structure::NodeIndexandEdge{0, ei2};

    //std::cout << ei1 << " " << ei2 << std::endl;

    msg_ = std::make_shared<motif_data_structure::MotifStateGraph>();
    msg_->add_state(m->get_state());
    if(parameters_.starting_helix != "") {
        auto h = rm_.motif_state(parameters_.starting_helix);
        msg_->add_state(h, start_.node_index, start_.edge_index);
        start_ = data_structure::NodeIndexandEdge{1, 1};
    }
}

void
DesignRNAScaffold::check_bp(
        String const & name,
        structure::RNAStructureOP struc,
        String const & type) {

    auto bps = structure::BasepairOPs();
    try { bps = struc->get_basepair(name); }
    catch(std::runtime_error) {
        LOG_ERROR << "cannot find " + type  + " basepair " + name; exit(0);
    }

    if(bps[0]->bp_type() != "cW-W" && parameters_.no_basepair_checks) {
        bps[0]->bp_type("cW-W");
        LOG_WARNING << "basepair " + name + " is not a watson and crick bp, but forcing it to be, this may be bad";
    }

    if(bps[0]->bp_type() != "cW-W") {
        LOG_ERROR << "basepair " + name + " is not a watson and crick bp "; exit(0);
    }
}

motif_search::SearchOP
DesignRNAScaffold::_setup_search() {
    auto search = motif_search::SearchOP(nullptr);

    if(parameters_.search_type == "path_finding") {
        LOGI << "Using Path_Finding search";
        using namespace motif_search::path_finding;
        auto scorer = std::make_shared<AstarScorer>();
        auto selector = default_selector();
        auto filter = std::make_shared<motif_search::NoExclusionFilter>();
        auto pf_search = std::make_shared<Search>(scorer, selector, filter);
        search = motif_search::SearchOP(pf_search->clone());
    }

    if(parameters_.search_type == "exhaustive") {
        LOGI << "Using Exhustive search";
        using namespace motif_search::exhaustive;
        auto score_factory = ScorerFactory();
        auto scorer = score_factory.get_scorer(parameters_.exhaustive_scorer);
        //TODO setup default sol topology if user does not specify motif_path
        auto sol_template = std::make_shared<motif_search::SolutionTopologyTemplate>();
        if(parameters_.motif_path != "") {
            sol_template = _setup_sol_template_from_path(parameters_.motif_path);
        }
        auto factory = motif_search::SolutionToplogyFactory();
        auto sol_toplogy = factory.generate_toplogy(*sol_template);
        auto filter = _setup_sol_filter(parameters_.solution_filter);
        auto e_search = std::make_shared<Search>(scorer, *sol_toplogy, filter);
        search = motif_search::SearchOP(e_search->clone());
    }

    if(parameters_.search_type == "mc") {
        LOGI << "Using Monte Carlo search";
        using namespace motif_search::monte_carlo;
        auto scorer = ScorerOP(nullptr);
        if(parameters_.mc_scorer == "default") {
            scorer = std::make_shared<DefaultScorer>();
        }
        else if(parameters_.mc_scorer == "scaled_scorer") {
            scorer = std::make_shared<ScaledScorer>(parameters_.scaled_score_d, parameters_.scaled_score_r);
        }

        if(parameters_.motif_path == "") {
            LOGE << "must include a motif path when using Monte Carlo search"; exit(0);
        }
        auto sol_template = _setup_sol_template_from_path(parameters_.motif_path);
        auto factory = motif_search::SolutionToplogyFactory();
        factory.set_option_value("max_helix_size", parameters_.max_helix_length);
        factory.set_option_value("min_helix_size", parameters_.min_helix_length);
        auto sol_toplogy = factory.generate_toplogy(*sol_template);
        auto filter = _setup_sol_filter("RemoveDuplicateHelices");
        auto e_search = std::make_shared<Search>(scorer, *sol_toplogy, filter);
        search = motif_search::SearchOP(e_search->clone());
    }

    // setup parameters
    search->set_option_value("accept_score", parameters_.search_cutoff);
    search->set_option_value("max_size", parameters_.search_max_size);


    return search;

}

motif_search::ProblemOP
DesignRNAScaffold::_setup_problem() {

    auto start_bp = msg_->get_node(start_.node_index)->data()->get_end_state(start_.edge_index);
    auto end_bp = msg_->get_node(end_.node_index)->data()->get_end_state(end_.edge_index);

    bool target_an_aligned_end = false;
    if(end_.edge_index == msg_->get_node(end_.node_index)->data()->block_end_add()) {
        target_an_aligned_end = true;
    }

    auto lookup = std::make_shared<util::StericLookupNew>();
    for(auto const & n : *msg_) {
        lookup->add_points(n->data()->cur_state->beads());
    }

    return std::make_shared<motif_search::Problem>(start_bp, end_bp, lookup, target_an_aligned_end);


}

motif_search::SolutionTopologyTemplateOP
DesignRNAScaffold::_setup_sol_template_from_path(
        String const & motif_path) {
    auto sol_template = std::make_shared<motif_search::SolutionTopologyTemplate>();
    auto spl = base::split_str_by_delimiter(motif_path, ",");
    int i = 0;
    auto lib_names = resources::MotifStateSqliteLibrary::get_libnames();
    for(auto const & e : spl) {
        if(lib_names.find(e) != lib_names.end()) {  // is a library
            if(i == 0) { sol_template->add_library(e); }
            else       { sol_template->add_library(e, data_structure::NodeIndexandEdge{i-1, 1}); }
        }

        else if(e == "tc_hairpin_hairpin") {
            std::cout << "made it" << std::endl;
        }

        else if(new_motif_ensembles_.find(e) != new_motif_ensembles_.end()) { // found user specified ensemble
            if(i == 0) { sol_template->add_ensemble(new_motif_ensembles_[e]); }
            else       { sol_template->add_ensemble(new_motif_ensembles_[e], data_structure::NodeIndexandEdge{i-1, 1}); }
        }
        else {
            // TODO should try and get motif from every end!
            auto ms = rm_.motif_state(e);
            if(i == 0) { sol_template->add_motif_state(ms); }
            else       { sol_template->add_motif_state(ms, data_structure::NodeIndexandEdge{i-1, 1}); }
        }
        i++;
    }

    return sol_template;
}

motif_search::SolutionFilterOP
DesignRNAScaffold::_setup_sol_filter(
        String const & solution_filter_name) {
    auto sf = motif_search::SolutionFilterOP(nullptr);
    if     (solution_filter_name == "NoFilter") {
        auto new_sf = std::make_shared<motif_search::NoExclusionFilter>();
        sf = motif_search::SolutionFilterOP(new_sf->clone());
    }
    else if (solution_filter_name == "RemoveDuplicateHelices") {
        auto new_sf = std::make_shared<motif_search::RemoveDuplicateHelices>();
        sf = motif_search::SolutionFilterOP(new_sf->clone());
    }
    else {
        LOG_ERROR << "not a valid solution filter name: " << solution_filter_name; exit(0);
    }

    return sf;
}

void
DesignRNAScaffold::_record_solution(
        motif_data_structure::MotifGraphOP mg,
        motif_search::SolutionOP design_sol,
        sequence_optimization::OptimizedSequenceOP sequence_opt_sol,
        int hit_count,
        int design_num,
        int sequence_opt_num) {
    score_out_ << design_num << "," << design_sol->score << "," << mg->designable_sequence() << ",";
    score_out_ << mg->dot_bracket() << "," << motif_names_ << ",";

    if(parameters_.dump_pdbs) {
        mg->to_pdb("design."+std::to_string(design_num)+".pdb", 1, 1);
    }

    if(sequence_opt_sol == nullptr) {
        score_out_ << ",,," << std::endl;
        out_ << mg->to_str() << std::endl;
        return;
    }
    else {
        score_out_ << sequence_opt_sol->sequence << "," << sequence_opt_sol->dist_score << ",";
        score_out_ << sequence_opt_sol->eterna_score << "," << hit_count << std::endl;
    }

    //score_out_ << "design_num,design_score,design_sequence,design_structure,motifs_uses,opt_num,";
    //score_out_ << "opt_sequence,opt_score,eterna_score" << std::endl;

}

void
DesignRNAScaffold::_fix_flex_helices_mtype(
        motif_data_structure::MotifGraphOP mg) {
    for(auto & n : *mg) {
        if(n->data()->name().substr(0, 5) == "HELIX") {
            n->data()->mtype(util::MotifType::HELIX);
        }
    }
}

String const &
DesignRNAScaffold::_get_motif_names(
        motif_data_structure::MotifGraphOP mg) {
    motif_names_ = "";
    for (auto const & n : *mg) {
        motif_names_ += n->data()->name() + ";";
    }
    return motif_names_;
}

void
DesignRNAScaffold::_build_new_ensembles(
        String const & new_ensemble_path) {
    auto lines = base::get_lines_from_file(new_ensemble_path);
    auto ensemble = motif::MotifStateEnsembleOP(nullptr);
    auto ms = motif::MotifStateOP(nullptr);
    auto all_ms = motif::MotifStateOPs();
    auto scores = Floats();
    for(auto const & l : lines) {
        auto spl = base::split_str_by_delimiter(l, " ");
        if(spl.size() != 3) { continue;}
        all_ms.resize(0); scores.resize(0);
        for(auto const & m_name : base::split_str_by_delimiter(spl[2], ",")) {
            ms = rm_.motif_state(m_name, "", "");
            for(auto const & end_name : ms->end_names()) {
                try {
                    all_ms.push_back(rm_.motif_state(m_name, "", end_name));
                    scores.push_back(0);
                }
                catch(...) { continue; }
            }
        }
        new_motif_ensembles_[spl[0]] = std::make_shared<motif::MotifStateEnsemble>(all_ms, scores);
    }

}

DesignRNAScaffold::EnsembleConversionResultsOP
DesignRNAScaffold::_get_mseg(
        motif_data_structure::MotifGraphOP mg) {

    auto index_hash = std::map<int, int>();
    auto mseg = std::make_shared<motif_data_structure::MotifStateEnsembleGraph>();
    int i = 0, j = 0;
    for(auto const & n : *mg) {
        // build / get motif state ensemble
        auto mse = motif::MotifStateEnsembleOP(nullptr);

        if(n->data()->mtype() == util::MotifType::HELIX) {
            // is not a basepair step
            if(n->data()->residues().size() > 4) {
                LOG_ERROR << "supplied a helix motif: "+n->data()->name() +" that is not a basepair step, this is no supported";
                exit(0);
            }

            try {
                mse = rm_.motif_state_ensemble(n->data()->end_ids()[0]);
            }
            catch(resources::ResourceManagerException const & e) {
                LOG_ERROR << "cannot find motif state ensemble for basepair with id: " + n->data()->end_ids()[0] <<
                          "check to make sure its a Watson-Crick basepair step";
                exit(0);
            }
        }
        else {
            mse = std::make_shared<motif::MotifStateEnsemble>(n->data()->get_state());
        }

        if(i == 0) {
            j = mseg->add_ensemble(*mse);
        }
        else {
            int pi = index_hash[mg->parent_index(n->index())];
            int pie = mg->parent_end_index(n->index());
            j = mseg->add_ensemble(*mse, data_structure::NodeIndexandEdge{pi, pie});
        }

        index_hash[n->index()] = j;
        i++;

    }
    return std::make_shared<EnsembleConversionResults>(mseg, index_hash);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int
main(
        int argc,
        const char ** argv) {

    //must add this for all apps!
    std::set_terminate(base::print_backtrace);

    //start logging
    base::init_logging();
    //load extra motifs being used
    String base_path = base::base_dir() + "/apps/simulate_tectos/resources/";
    resources::Manager::instance().add_motif(base_path + "GAAA_tetraloop");

    auto app = DesignRNAScaffold();
    app.setup_options();
    CLI11_PARSE(app.app_, argc, argv);
    app.run();

    return 0;

}
