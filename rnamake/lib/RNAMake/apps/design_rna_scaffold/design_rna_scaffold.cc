//
// Created by Joseph Yesselman on 3/9/19.
//

#include "base/backtrace.hpp"
#include "base/log.h"
#include "design_rna_scaffold/design_rna_scaffold.h"

#include <motif_search/solution_topology.h>
#include <motif_search/path_finding/search.h>
#include <motif_search/exhaustive/search.h>
#include <motif_search/monte_carlo/search.h>


DesignRNAScaffold::DesignRNAScaffold():
        base::Application(),
        rm_(resources::Manager::instance()),
        parameters_(Parameters()) {

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// app functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
DesignRNAScaffold::setup_options() {
    // core inputs
    // from pdb
    add_option("pdb", String(""), base::OptionType::STRING, false);
    add_option("start_bp", String(""), base::OptionType::STRING, false);
    add_option("end_bp", String(""), base::OptionType::STRING, false);
    // from motif graph (not used much)
    add_option("mg", String(""), base::OptionType::STRING, false);

    // common options
    add_option("designs", 1, base::OptionType::INT, false);
    add_option("dump_pdbs", false, base::OptionType::BOOL, false);
    add_option("dump_scaffold_pdbs", false, base::OptionType::BOOL, false);
    add_option("helix_type", "ideal_helices", base::OptionType::STRING, false);
    add_option("log_level", "info", base::OptionType::STRING, false);
    add_option("search_type", "path_finding", base::OptionType::STRING, false);
    add_option("search_cutoff", 5.0f, base::OptionType::FLOAT, false);
    add_option("search_max_size", 999999, base::OptionType::INT, false);
    add_option("skip_sequence_optimization", false, base::OptionType::BOOL, false);
    add_option("sequence_opt_cutoff", 5.0f, base::OptionType::FLOAT, false);
    add_option("solution_path", "", base::OptionType::STRING, false);
    add_option("out_file", "default.out", base::OptionType::STRING, false);
    add_option("score_file", "default.scores", base::OptionType::STRING, false);
    add_option("solution_filter", "NoFilter", base::OptionType::STRING, false);

    // less common options
    add_option("no_basepair_checks", false, base::OptionType::BOOL, false);
    add_option("no_out_file", false, base::OptionType::BOOL, false);
    add_option("max_helix_length", 10, base::OptionType::INT, false);
    add_option("min_helix_length", 3, base::OptionType::INT, false);
    add_option("starting_helix", "", base::OptionType::STRING, false);
    add_option("ending_helix", "", base::OptionType::STRING, false);
    add_option("all_designs", false, base::OptionType::BOOL, false);
    add_option("flip_start_bp", false, base::OptionType::BOOL, false);
    add_option("flip_end_bp", false, base::OptionType::BOOL, false);
    add_option("motif_path", "", base::OptionType::STRING, false);
    add_option("new_ensembles", "", base::OptionType::STRING, false);

    // scoring related options
    add_option("exhaustive_scorer", "default", base::OptionType::STRING, false);
    add_option("mc_scorer", "default", base::OptionType::STRING, false);
    add_option("scaled_score_d", 1.0f, base::OptionType::FLOAT, false);
    add_option("scaled_score_r", 2.0f, base::OptionType::FLOAT, false);

}

void
DesignRNAScaffold::parse_command_line(
        int argc,
        const char ** argv) {
    base::Application::parse_command_line(argc, argv);

    // core inputs
    parameters_.pdb       = get_string_option("pdb");     parameters_.start_bp   = get_string_option("start_bp");
    parameters_.end_bp    = get_string_option("end_bp");  parameters_.mg         = get_string_option("mg");
    // common options
    parameters_.designs   = get_int_option("designs");
    parameters_.dump_pdbs = get_bool_option("dump_pdbs"); parameters_.dump_scaffold_pdbs = get_bool_option("dump_scaffold_pdbs");
    parameters_.out_file = get_string_option("out_file"); parameters_.score_file = get_string_option("score_file");
    parameters_.solution_filter = get_string_option("solution_filter");
    parameters_.search_type = get_string_option("search_type");
    parameters_.motif_path = get_string_option("motif_path");
    parameters_.all_designs = get_bool_option("all_designs");
    parameters_.skip_sequence_optimization = get_bool_option("skip_sequence_optimization");
    // less common options
    parameters_.no_basepair_checks         = get_bool_option("no_basepair_checks");
    parameters_.starting_helix = get_string_option("starting_helix");
    parameters_.ending_helix = get_string_option("ending_helix");
    parameters_.search_cutoff = get_float_option("search_cutoff");
    parameters_.search_max_size = get_int_option("search_max_size");
    parameters_.new_ensembles = get_string_option("new_ensembles");
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
    score_out_ << "opt_sequence,opt_score,eterna_score" << std::endl;

    search_ = _setup_search();
    problem_ = _setup_problem();

    auto mg = msg_->to_motif_graph();
    auto sol_mg = motif_data_structure::MotifGraphOP(nullptr);
    mg->increase_level();

    auto sol = motif_search::SolutionOP(nullptr);
    search_->setup(problem_);
    int i = 0;
    while(!search_->finished()) {
        sol = search_->next();
        std::cout << sol << std::endl;
        if(sol == nullptr) { break; }

        sol_mg = sol->graph->to_motif_graph();
        if(sol_mg->size() == 0) { break; }
        _get_motif_names(sol_mg);

        if(i == 0) {
            LOG_INFO << "found a solution: " << motif_names_;
        }
        else if(i % 10 == 0) {
            LOG_INFO << "found " << i << " solutions ";
        }

        mg->add_motif_graph(*sol_mg, start_.node_index, start_.edge_index);
        sol_mg->to_pdb("test.pdb", 1, 1);
        auto mg_copy = std::make_shared<motif_data_structure::MotifGraph>(*mg);
        mg_copy->add_connection(mg->last_node()->index(), end_.node_index,
                                mg->last_node()->data()->end_name(1),
                                mg->get_node(end_.node_index)->data()->end_name(end_.edge_index));
        _fix_flex_helices_mtype(mg_copy);
        mg_copy->replace_ideal_helices();

        _record_solution(mg_copy, sol, nullptr, i, 0);
        mg->remove_level(1);

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
    m->to_pdb("start.pdb");

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
    String base_path = base::base_dir() + "/rnamake/lib/RNAMake/apps/simulate_tectos/resources/";
    resources::Manager::instance().add_motif(base_path + "GAAA_tetraloop");

    auto app = DesignRNAScaffold();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;

}




















