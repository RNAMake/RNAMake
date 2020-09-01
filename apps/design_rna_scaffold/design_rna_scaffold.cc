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

DesignRNAScaffold::DesignRNAScaffold():
        base::Application(),
        rm_(resources::Manager::instance()),
        app_("DesignRNAScaffold") {


}

String
valid_pdb(String& path) {
    auto ending = path.substr(path.size()-4);
    return ending == ".pdb" ? String{""} : String{"the file specified by --pdb must end in .pdb"};
}

String
valid_bp(String& bp) {
    auto bp_pattern = std::regex("[A-Z][0-9]{*}-[A-Z][0-9]{*}");
    std::regex_match(bp.begin(),bp.end(),bp_pattern);
    std::cout<<bp_pattern.mark_count()<<std::endl;
    return String{"bad"};
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// app functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
DesignRNAScaffold::setup_options() {

    ////////////////////////////////////////////////////////////////////////////////
    // Core Inputs
    ////////////////////////////////////////////////////////////////////////////////

    app_.add_option_group("Core Inputs");
    app_.add_option_group("I/O Options");
    app_.add_option_group("Search Parameters");
    app_.add_option_group("Scoring Paramters");
    app_.add_option_group("Sequence Optimization Parameters");
    app_.add_option_group("TBD");

    // Core Inputs
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    app_.add_option("--pdb",
                    parameters_.core.pdb,
                    "path to a PDB file with input RNA structure")
                    ->required()
                    ->check(CLI::ExistingFile&CLI::Validator(valid_pdb,"ends in .pdb","valid_pdb"))
                    ->group("Core Inputs");
<<<<<<< HEAD

    app_.add_option("--start_bp",
                        parameters_.core.start_bp,
                        "starting basepair to be used in structure format: [CHAIN ID][NT1 NUM] - [CHAIN ID][NT2 NUM]")
=======

    app_.add_option("--end_bp",
                    parameters_.end_bp,
                    "ending basepair to be used in structure format: [CHAIN ID][NT1 NUM]-[CHAIN ID][NT2 NUM]")
>>>>>>> devel
                    ->required()
                    ->check(CLI::Validator(valid_bp,"format [CHAIN ID][NT1 NUM]-[CHAIN ID][NT2 NUM]","valid_bp"))
                    ->group("Core Inputs");

<<<<<<< HEAD
    app_.add_option("--end_bp",
                    parameters_.core.end_bp,
                    "ending basepair to be used in structure format: [CHAIN ID][NT1 NUM] - [CHAIN ID][NT2 NUM]")
=======
    app_.add_option("--start_bp",
                    parameters_.start_bp,
                    "starting basepair to be used in structure format: [CHAIN ID][NT1 NUM]-[CHAIN ID][NT2 NUM]")
>>>>>>> devel
                    ->required()
                    ->group("Core Inputs");

    app_.add_option("--designs",
                        parameters_.core.designs,
                        "number of designs to create. Default is 1")
                    ->default_val(1)
                    ->check(CLI::PositiveNumber)
                    ->group("Core Inputs");
    //CLI::Validator(std::function<std::string(std::string&)>,);
    ////////////////////////////////////////////////////////////////////////////////
    // File Options
    ////////////////////////////////////////////////////////////////////////////////

    app_.add_option_group("File Options");

    // I/O Options
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    app_.add_option("--out_file",
                    parameters_.io.out_file,
                    "output file that contains all information to rebuild solutions")
                    ->default_val("default.out")
                    ->group("I/O Options");

    app_.add_option("--score_file",
                    parameters_.io.score_file,
                    parameters_.io.score_file,
                    "name of output file containining scoring information for design")
                    ->default_val("default.scores")
                    ->group("File Options");

    app_.add_flag("--dump_pdbs",
                    parameters_.io.dump_pdbs,
                    "flag to dump intermediate pdbs TODO")
                    ->group("I/O Options");

    app_.add_flag("--dump_scaffold_pdbs",
<<<<<<< HEAD
                    parameters_.io.dump_scaffold_pdbs,
                    "flag to output pdbs of just the design scaffold WITHOUT initial RNA structure useful for really big structures")
                    ->group("I/O Options");

    app_.add_option("--new_ensembles",
                    parameters_.io.new_ensembles_file,
                    "flag to include new structural ensembles")
                    ->default_val("")
                    ->group("I/O Options");

    app_.add_flag("--no_out_file",
                  parameters_.io.no_out_file,
                  "if you only want the summary and not the actual structures")
                  ->default_val(false)
                  ->group("I/O Options");

    // Search Options
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

=======
                    parameters_.dump_scaffold_pdbs,
                    "flag to output intermediate scaffold pdbs in the current directory")
                    ->group("File Options");
    app_.add_option("--score_file",
                    parameters_.score_file,
                    "name of output file containining scoring information for design")
                    ->default_val("default.scores")
                    ->group("File Options");
    ////////////////////////////////////////////////////////////////////////////////
    // Search Parameters
    ////////////////////////////////////////////////////////////////////////////////
    app_.add_option_group("Search Parameters");

>>>>>>> devel
    app_.add_option("--search_type",
                    parameters_.search.type,
                    "search type for traversing motif space")
                    ->default_val("path_finding")
                    ->check(CLI::IsMember(std::set<String>{"path_finding", "exhaustive", "mc"}))
                    ->group("Search Parameters");

    app_.add_option("--search_cutoff",
                    parameters_.search.cutoff,
                    "TODO")
                    ->default_val(5.0f)
                    ->group("Search Parameters");

    app_.add_option("--search_max_size",
                    parameters_.search.max_size,
                    "maximum number of steps for a design search")
                    ->default_val(999999)
                    ->check(CLI::PositiveNumber) //TODO put a limit to this?? CJ 08/20
                    ->group("Search Parameters");

<<<<<<< HEAD
=======
    app_.add_flag("--skip_sequence_optimization",
                    parameters_.skip_sequence_optimization,
                    "flag to skip sequence optimization of the design")
                    ->group("Search Parameters");

    app_.add_option("--solution_filter",
                    parameters_.solution_filter,
                    "TODO")
                    ->default_val("NotFilter")
                    ->check(CLI::IsMember(std::set<String>{"NoFilter","RemoveDuplicateHelices"}))
                    ->group("Search Parameters");

>>>>>>> devel
    app_.add_flag("--no_basepair_checks",
                  parameters_.search.no_basepair_checks,
                  "flag to disable basepair checks")
                  ->group("Search Parameters");

    app_.add_option("--max_helix_length",
                    parameters_.search.max_helix_length,
                    "maximum number of basepairs in a solution helix")
                    ->default_val(99)
                    ->group("Search Parameters");

    app_.add_option("--min_helix_length",
                    parameters_.search.min_helix_length,
                    "minimum number of basepairs in a solution helix")
                    ->default_val(4)
                    ->group("Search Parameters");
    app_.add_option("--starting_helix",
                    parameters_.search.starting_helix,
                    "starting helix for design solution. Format = [TODO]")
                    ->default_val("")
                    ->group("Search Parameters");

    app_.add_option("--ending_helix",
                    parameters_.search.ending_helix,
                    "ending helix for design solution. Format = [TODO]")
                    ->default_val("")
                    ->group("Search Parameters");
    ////////////////////////////////////////////////////////////////////////////////
    // Scoring Paramters
    ////////////////////////////////////////////////////////////////////////////////
    app_.add_option_group("Scoring Paramteres");

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
    ////////////////////////////////////////////////////////////////////////////////
    // TBD
    ////////////////////////////////////////////////////////////////////////////////
    app_.add_option_group("TBD");

    app_.add_flag("--skip_thermo_fluc",
                  parameters_.skip_thermo_fluc,
                  "flag to skip thermo fluc calculate to evaluate stability")
                  ->group("TBD");

    app_.add_option("--motif_path",
                    parameters_.search.motif_path,
                    "TBD")
                    ->default_val("")
                    ->group("Search Parameters");

    app_.add_option("--solution_filter",
                    parameters_.search.solution_filter,
                    "TODO")
                    ->default_val("NotFilter")
                    ->check(CLI::IsMember(std::set<String>{"NoFilter","RemoveDuplicateHelices"}))
                    ->group("Search Parameters");

<<<<<<< HEAD
    app_.add_option("--exhaustive_scorer",
                    parameters_.search.exhaustive_scorer,
                    "TODO")
                    ->default_val("default")
                    ->group("Scoring Parameters");

    app_.add_option("--mc_scorer",
                    parameters_.search.mc_scorer,
                    "TODO")
                    ->default_val("default")
                    ->check(CLI::IsMember(std::set<String>{"default","scaled_scorer"}))
                    ->group("Scoring Parameters");

    app_.add_option("--scaled_score_d",
                    parameters_.search.scaled_score_d,
                    "TODO")
                    ->default_val(1.0f)
                    ->group("Scoring Parameters");

    app_.add_option("--scaled_score_r",
                    parameters_.search.scaled_score_r,
                    "TODO")
                    ->default_val(2.0f)
                    ->group("Scoring Parameters");

    // Sequence Optimization Options
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    app_.add_flag("--skip_sequence_optimization",
                    parameters_.seq_opt.skip,
                    "flag to skip sequence optimization of the design")
                    ->group("Search Parameters");

    app_.add_flag("--skip_thermo_fluc",
                  parameters_.thermo_fluc.skip,
                  "flag to skip thermo fluc calculate to evaluate stability")
                  ->group("TBD");


    // less common options





=======

>>>>>>> devel
}

void
DesignRNAScaffold::parse_command_line(
        int argc,
        const char ** argv) {

 }

void
DesignRNAScaffold::run() {
    // sets up all variables required
    setup();

    LOG_INFO << "###################";
    LOG_INFO << "# starting search #";
    LOG_INFO << "###################";
    while(true) {
        auto found = _get_motif_graph_solution();
        if(!found) {
            break;
        }
        if(parameters_.seq_opt.skip) {
            design_num_ += 1;
            _record_solution(mg_w_sol_, solution_, nullptr, 0, design_num_, 0);
            continue;
        }
        mg_w_sol_->replace_ideal_helices();
        // remap node indexes with new base pair step motifs
        auto last_m = mg_w_sol_->get_node(end_.node_index)->connections()[end_.edge_index]->partner(end_.node_index);
        auto new_start = data_structure::NodeIndexandEdge{last_m->index(), 1};

        auto opt_seq_scorer = std::make_shared<sequence_optimization::InternalTargetScorer>(
            new_start.node_index, new_start.edge_index, end_.node_index, end_.edge_index,
            problem_->target_an_aligned_end);
        auto sols = seq_optimizer_->get_optimized_sequences(mg_w_sol_, opt_seq_scorer);
        if(sols.empty()) {
            continue;
        }
        mg_w_sol_->replace_helical_sequence(sols[0]->sequence);
        _record_solution(mg_w_sol_, solution_, sols[0], 0, design_num_, 0);
        if(parameters_.core.designs <= design_num_) {
            LOG_INFO << "found " << design_num_ << " designs, if you want more please use --designs num_of_design";
            exit(0);
        }
    }

    exit(0);

    /*while(!search_->finished()) {
        try {

            if(parameters_.skip_sequence_optimization) {
                i += 1;
                _record_solution(mg_copy, sol, nullptr, 0, i, 0);
                mg->remove_level(1);
                continue;
            }
            mg_copy->replace_ideal_helices();
            auto last_m = mg_copy->get_node(end_.node_index)->connections()[end_.edge_index]->partner(end_.node_index);
            auto new_start = data_structure::NodeIndexandEdge{last_m->index(), 1};

            for (int j = 0; j < 1; j++) {

                // sequence optimization
                auto opt_seq_scorer = std::make_shared<sequence_optimization::InternalTargetScorer>(
                        new_start.node_index, new_start.edge_index, end_.node_index, end_.edge_index,
                        problem_->target_an_aligned_end);
                auto sols = seq_optimizer->get_optimized_sequences(mg_copy, opt_seq_scorer);
                //std::cout << seq_optimizer->get_float_option("cutoff") << std::endl;
                //std::cout << sols[0]->dist_score << std::endl;
                mg_copy->replace_helical_sequence(sols[0]->sequence);


                // thermo sim
                auto r = _get_mseg(mg_copy);

                auto start = data_structure::NodeIndexandEdge{r->index_hash[new_start.node_index],
                                                              new_start.edge_index};
                auto end = data_structure::NodeIndexandEdge{r->index_hash[end_.node_index], end_.edge_index};

                sim->setup(*r->mseg, start, end);
                auto count = 0;
                for (int s = 0; s < 1000000; s++) {
                    if(sim->next()) {
                        count += 1;
                    }

                }
                _record_solution(mg_copy, sol, sols[0], count, i, j);

            }


            mg->remove_level(1);
        }
        catch(std::runtime_error const & e) {
            LOG_WARNING << "error was caught:" << e.what();
            mg->remove_level(1);
            continue;
        }

        i++;
        if(parameters_.designs <= i) {
            LOG_INFO << "found " << i << " designs, if you want more please use --designs num_of_design";
            exit(0);
        }
    }

    if(i == 0) {
        LOG_INFO << "no solutions found";
    }
    else {
        LOG_INFO << "no more solutions found";
        LOG_INFO << "found a total of " << i << " solution(s)";
    }*/

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
DesignRNAScaffold::setup() {
    LOG_INFO << "#########";
    LOG_INFO << "# setup #";
    LOG_INFO << "#########";

    if(!parameters_.io.new_ensembles_file.empty()) {
        _build_new_ensembles(parameters_.io.new_ensembles_file);
    }
    // setups up start_, end_ and msg_;
    _setup_from_pdb();

    mg_ = msg_->to_motif_graph();
    mg_->set_option_value("sterics", false);
    mg_->increase_level();

    // sets up outputing files
    LOG_INFO << "out file that contains all build solutions -> " + parameters_.io.out_file + " can be set with --out_file";
    LOG_INFO << "score file that contains summary information about solutions -> " + parameters_.io.score_file + " can be set with --score_file";
    out_.open(parameters_.io.out_file);
    score_out_.open(parameters_.io.score_file);

    score_out_ << "design_num,design_score,design_sequence,design_structure,motifs_uses,opt_num,";
    score_out_ << "opt_sequence,opt_score,eterna_score,hit_count" << std::endl;

    // setup search
    search_ = _setup_search();
    problem_ = _setup_problem();
    search_->setup(problem_);

    //sequence optimziation setup
    seq_optimizer_ = std::make_shared<sequence_optimization::SequenceOptimizer3D>();
    seq_optimizer_->set_option_value("steps", 1000);

    if(! parameters_.seq_opt.skip) {
        LOG_INFO << "sequence optimization cutoff -> " << seq_optimizer_->get_float_option("cutoff");
        LOG_INFO << "sequence optimization steps -> " << seq_optimizer_->get_int_option("steps");
    }
    else {
        LOG_INFO << "set --skip_sequence_optimization -> skipping sequence optimization step!";
        LOG_WARNING << "Note: skipping sequence_optimization may result in more chain breaks in outputed pdbs";
    }

    //thermo sim setup
    auto thermo_scorer = std::make_shared<thermo_fluctuation::graph::FrameScorer>();
    auto sterics = std::make_shared<thermo_fluctuation::graph::sterics::NoSterics>();
    thermo_sim_ = std::make_shared<thermo_fluctuation::graph::Simulation>(thermo_scorer, sterics);

}

bool
DesignRNAScaffold::_get_motif_graph_solution() {
    solution_ = search_->next();
    // search could not find a solution
    if (solution_ == nullptr) {
        return false;
    }
    auto sol_mg = solution_->graph->to_motif_graph();
    // for some reason solution is size 0 this reall shouldnt happen
    if (sol_mg->size() == 0) {
        return false;
    }
    _get_motif_names(sol_mg);

    if (design_num_ == 0) {
        LOG_INFO << "found a solution: " << motif_names_ << " with score: " << solution_->score;
    } else if (design_num_ % 10 == 0) {
        LOG_INFO << "found " << design_num_ << " solutions ";
    }

    mg_->add_motif_graph(*sol_mg, start_.node_index, start_.edge_index);
    mg_w_sol_ = std::make_shared<motif_data_structure::MotifGraph>(*mg_);
    mg_w_sol_->add_connection(mg_->last_node()->index(), end_.node_index,
                              mg_->last_node()->data()->end_name(1),
                              mg_->get_node(end_.node_index)->data()->end_name(end_.edge_index));
    _fix_flex_helices_mtype(mg_w_sol_);
    mg_->remove_level(1);
    return true;
}

void
DesignRNAScaffold::_setup_from_pdb() {
    auto struc = rm_.get_structure(parameters_.core.pdb, "scaffold");
    LOG_INFO << "loaded pdb from file: " << parameters_.core.pdb;
    if(struc->ends().empty()) {
        LOG_ERROR << "pdb has no ends available to build off of";
        exit(0);
    }
    auto end_names = String();
    for(auto const & end : struc->ends()) {
          end_names += end->name() + " ";
    }
    LOG_INFO << "available ends to build off: " << end_names;
    check_bp(parameters_.core.start_bp, struc, "start");
    check_bp(parameters_.core.end_bp,  struc, "end");

    // TODO allow for construction of motif from RNA structure to allow for changes in base pair types
    rm_.add_motif(parameters_.core.pdb, "scaffold", util::MotifType::TCONTACT);
    auto m = rm_.motif("scaffold", "", parameters_.core.end_bp);
    int ei1, ei2;

    try {
        ei1 = m->get_end_index(parameters_.core.start_bp);
    }
    catch(structure::RNAStructureException const & e) {
        LOG_ERROR << "starting basepair: " << parameters_.core.start_bp << " was NOT FOUND!";
        exit(0);
    }
    LOG_INFO << "start basepair: " << parameters_.core.start_bp << " is found!";
    try {
        ei2 = m->get_end_index(parameters_.core.end_bp);
    }
    catch(structure::RNAStructureException const & e) {
        LOG_ERROR << "end basepair: " << parameters_.core.end_bp << " was NOT FOUND!";
        exit(0);
    }
    LOG_INFO << "end basepair: " << parameters_.core.start_bp << " is found!";

    start_ = data_structure::NodeIndexandEdge{0, ei1};
    end_   = data_structure::NodeIndexandEdge{0, ei2};

    msg_ = std::make_shared<motif_data_structure::MotifStateGraph>();
    msg_->add_state(m->get_state());
    if(! parameters_.search.starting_helix.empty()) {
        auto h = rm_.motif_state(parameters_.search.starting_helix);
        msg_->add_state(h, start_.node_index, start_.edge_index);
        start_ = data_structure::NodeIndexandEdge{1, 1};
    }
}

void
DesignRNAScaffold::check_bp(
        String const & name,
        structure::RNAStructureOP const & struc,
        String const & type) {

    auto bps = structure::BasepairOPs();
    try { bps = struc->get_basepair(name); }
    catch(std::runtime_error const & e) {
        LOG_ERROR << "cannot find " + type  + " basepair " + name;
        exit(0);
    }

    if(bps[0]->bp_type() != "cW-W" && parameters_.search.no_basepair_checks) {
        bps[0]->bp_type("cW-W");
        LOG_WARNING << "basepair " + name + " is not a watson and crick bp, but forcing it to be, this may be bad";
    }

    if(bps[0]->bp_type() != "cW-W") {
        LOG_ERROR << "basepair " + name + " is not a watson and crick bp ";
        exit(0);
    }

}

motif_search::SearchOP
DesignRNAScaffold::_setup_search() {
    auto search = motif_search::SearchOP(nullptr);

    LOGI << "search_type -> " + parameters_.search.type;
    if(parameters_.search.type == "path_finding") {
        using namespace motif_search::path_finding;
        auto scorer = std::make_shared<AstarScorer>();
        auto selector = default_selector();
        auto filter = std::make_shared<motif_search::NoExclusionFilter>();
        auto pf_search = std::make_shared<Search>(scorer, selector, filter);
        search = motif_search::SearchOP(pf_search->clone());
    }

    if(parameters_.search.type == "exhaustive") {
        using namespace motif_search::exhaustive;
        auto score_factory = ScorerFactory();
        auto scorer = score_factory.get_scorer(parameters_.search.exhaustive_scorer);
        //TODO setup default sol topology if user does not specify motif_path
        auto sol_template = std::make_shared<motif_search::SolutionTopologyTemplate>();
        if(parameters_.search.motif_path != "") {
            sol_template = _setup_sol_template_from_path(parameters_.search.motif_path);
        }
        auto factory = motif_search::SolutionToplogyFactory();
        auto sol_toplogy = factory.generate_toplogy(*sol_template);
        auto filter = _setup_sol_filter(parameters_.search.solution_filter);
        auto e_search = std::make_shared<Search>(scorer, *sol_toplogy, filter);
        search = motif_search::SearchOP(e_search->clone());
    }

    if(parameters_.search.type == "mc") {
        using namespace motif_search::monte_carlo;
        auto scorer = ScorerOP(nullptr);
        if(parameters_.search.mc_scorer == "default") {
            scorer = std::make_shared<DefaultScorer>();
        }
        else if(parameters_.search.mc_scorer == "scaled_scorer") {
            scorer = std::make_shared<ScaledScorer>(
                parameters_.search.scaled_score_d,
                parameters_.search.scaled_score_r);
        }
        if(parameters_.search.motif_path == "") {
            LOGE << "must include a motif path when using Monte Carlo search"; exit(0);
        }
        auto sol_template = _setup_sol_template_from_path(parameters_.search.motif_path);
        auto factory = motif_search::SolutionToplogyFactory();
        factory.set_option_value("max_helix_size", parameters_.search.max_helix_length);
        factory.set_option_value("min_helix_size", parameters_.search.min_helix_length);
        auto sol_toplogy = factory.generate_toplogy(*sol_template);
        auto filter = _setup_sol_filter("RemoveDuplicateHelices");
        auto e_search = std::make_shared<Search>(scorer, *sol_toplogy, filter);
        search = motif_search::SearchOP(e_search->clone());
    }

    // setup parameters
    search->set_option_value("accept_score", parameters_.search.cutoff);
    search->set_option_value("max_size", parameters_.search.max_size);

    LOGI << "search cutoff -> " + std::to_string(parameters_.search.cutoff);
    LOGI << "search max size (how many residues) -> " + std::to_string(parameters_.search.max_size);
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

    if(parameters_.io.dump_pdbs) {
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
