//
// Created by Joseph Yesselman on 3/9/19.
//

#include <CLI/CLI.hpp>

#include "base/backtrace.h"
#include "base/log.h"
#include "design_rna_scaffold/design_rna_scaffold.h"

#include <motif_search/solution_topology.h>
#include <motif_search/path_finding/search.h>
#include <motif_search/exhaustive/search.h>
#include <motif_search/monte_carlo/search.h>

#include <sequence_optimization/sequence_optimizer.h>

DesignRNAScaffold::DesignRNAScaffold () :
    base::Application(),
    rm_(resources::Manager::instance()),
    app_("DesignRNAScaffold") {}


// app functions  //////////////////////////////////////////////////////////////////////////////////

String
valid_pdb (String &path) {
    auto ending = path.substr(path.size() - 4);
    return ending == ".pdb" ? String{""} : String{"the file specified by --pdb must end in .pdb"};
}

String
valid_bp (String &bp) {
    const auto bp_pattern = std::regex("\\b[A-Z][0-9]*-[A-Z][0-9]*\\b");
    auto sm = std::smatch{};
    std::regex_match(bp, sm, bp_pattern);

    return sm.size() == 1 ? String{""} : String{bp + " is an invalid bp format"};
}

void
DesignRNAScaffold::setup_options () {

    app_.add_option_group("Core Inputs");
    app_.add_option_group("I/O Options");
    app_.add_option_group("Search Parameters");
    app_.add_option_group("Scoring Paramters");
    app_.add_option_group("Sequence Optimization Parameters");
    app_.add_option_group("Thermo Fluc Parameters");
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Core Inputs
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    app_.add_option("--pdb", parameters_.core.pdb, "path to a PDB file with input RNA structure")
        ->check(CLI::ExistingFile & CLI::Validator(valid_pdb, "ends in .pdb", "valid_pdb"))
        ->group("Core Inputs");

    app_.add_option("--start_bp",
                    parameters_.core.start_bp,
                    "starting basepair to be used in structure format: [CHAIN ID][NT1 NUM]-[CHAIN ID][NT2 NUM]")
        ->check(CLI::Validator(valid_bp,
                               "format [CHAIN ID][NT1 NUM]-[CHAIN ID][NT2 NUM]",
                               "valid_bp"))
        ->group("Core Inputs");

    app_.add_option("--end_bp",
                    parameters_.core.end_bp,
                    "ending basepair to be used in structure format: [CHAIN ID][NT1 NUM]-[CHAIN ID][NT2 NUM]")
        ->check(CLI::Validator(valid_bp,
                               "format [CHAIN ID][NT1 NUM]-[CHAIN ID][NT2 NUM]",
                               "valid_bp"))
        ->group("Core Inputs");

    //app_.add_option("--mg",parameters_.core.mg,"path to a motif graph file")
    //                ->check(CLI::ExistingFile)
    //                ->group("Core Inputs");

    app_.add_option("--designs",
                    parameters_.core.designs,
                    "number of designs to create. Default is 1")
        ->default_val(1)
        ->check(CLI::PositiveNumber)
        ->group("Core Inputs");

    app_.add_option("--log_level", parameters_.core.log_level, "level for global logging")
        ->check(CLI::IsMember(std::set<String>{"debug", "error", "fatal", "info", "verbose",
                                               "warn"}))
        ->default_val("info")
        ->group("Core Inputs");

    app_.add_option("--extra_pdbs", parameters_.core.extra_pdbs, ", deliminted list of other pdbs used in building")
        ->default_val("")
        ->group("Core Inputs");
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // I/O Options
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    app_.add_flag("--dump_intermediate_pdbs",
                  parameters_.io.dump_intermediate_pdbs,
                  "flag to dump intermediate pdbs")
        ->group("I/O Options");

    app_.add_flag("--dump_pdbs", parameters_.io.dump_pdbs, "TODO")
        ->group("I/O Options");

    app_.add_flag("--dump_scaffold_pdbs",
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

    app_.add_option("--out_file",
                    parameters_.io.out_file,
                    "output file that contains all information to rebuild solutions")
        ->default_val("default.out")
        ->group("I/O Options");

    app_.add_option("--score_file",
                    parameters_.io.score_file,
                    "name of output file containining scoring information for design")
        ->default_val("default.scores")
        ->group("I/O Options");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Search Parameters
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    app_.add_option("--ending_helix",
                    parameters_.search.ending_helix,
                    "ending helix for design solution. Format = [TODO]")
        ->default_val("")
        ->group("Search Parameters");

    app_.add_option("--exhaustive_scorer", parameters_.search.exhaustive_scorer, "TODO")
        ->default_val("default")
        ->group("Search Parameters");

    app_.add_option("--max_helix_length",
                    parameters_.search.max_helix_length,
                    "maximum number of basepairs in a solution helix")
        ->default_val(99)
        ->group("Search Parameters");

    app_.add_option("--mc_scorer", parameters_.search.mc_scorer, "TODO")
        ->default_val("default")
        ->check(CLI::IsMember(std::set<String>{"default", "scaled_scorer"}))
        ->group("Search Parameters");

    app_.add_option("--min_helix_length",
                    parameters_.search.min_helix_length,
                    "minimum number of basepairs in a solution helix")
        ->default_val(4)
        ->group("Search Parameters");

    app_.add_option("--motif_path", parameters_.search.motif_path, "TBD")
        ->group("Search Parameters");

    app_.add_flag("--no_basepair_checks",
                  parameters_.search.no_basepair_checks,
                  "flag to disable basepair checks")
        ->group("Search Parameters");

    app_.add_option("--scaled_score_d", parameters_.search.scaled_score_d, "TODO")
        ->default_val(1.0f)
        ->group("Search Parameters");

    app_.add_option("--scaled_score_r", parameters_.search.scaled_score_r, "TODO")
        ->default_val(2.0f)
        ->group("Search Parameters");

    app_.add_option("--search_cutoff", parameters_.search.cutoff, "TODO")
        ->default_val(7.5f)
        ->group("Search Parameters");

    app_.add_option("--search_max_size",
                    parameters_.search.max_size,
                    "maximum number of steps for a design search")
        ->default_val(999999)
        ->check(CLI::PositiveNumber) //TODO put a limit to this?? CJ 08/20
        ->group("Search Parameters");

    app_.add_option("--search_type",
                    parameters_.search.type,
                    "search type for traversing motif space")
        ->default_val("path_finding")
        ->check(CLI::IsMember(std::set<String>{"path_finding", "exhaustive", "mc"}))
        ->group("Search Parameters");

    app_.add_option("--solution_filter", parameters_.search.solution_filter, "TODO")
        ->default_val("RemoveDuplicateHelices")
        ->check(CLI::IsMember(std::set<String>{"NoFilter", "RemoveDuplicateHelices"}))
        ->group("Search Parameters");

    app_.add_option("--starting_helix",
                    parameters_.search.starting_helix,
                    "starting helix for design solution. Format = [TODO]")
        ->default_val("")
        ->group("Search Parameters");

    app_.add_flag("--no_sterics",
                  parameters_.search.no_sterics,
                  "turns off sterics checks againsts supplied RNA structure")
        ->group("Search Parameters");

    app_.add_flag("--only_tether_opt",
                  parameters_.search.only_tether_opt,
                  "ignore supplied structure other than sterics")
        ->group("Search Parameters");

    app_.add_option("--search_max_motifs", parameters_.search.max_motifs, "TODO")
        ->default_val(999)
        ->group("Search Parameters");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sequence Optimization Options
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    app_.add_flag("--skip_sequence_optimization",
                  parameters_.seq_opt.skip,
                  "flag to skip sequence optimization of the design")
        ->group("Sequence Optimization Parameters");

    app_.add_option("--sequences_per_design",
                    parameters_.seq_opt.sequences_per_design,
                    "number of sequences to try per motif design")
        ->default_val(1)
        ->group("Sequence Optimization Parameters");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Thermo fluc Options
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    app_.add_flag("--thermo_fluc",
                  parameters_.thermo_fluc.perform,
                  "run thermo fluc procedure to estimate thermo fluc of helices")
        ->group("Thermo Fluc Parameters");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Setting some global parameters/variables for the parser
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //app_.get_option("--pdb");//->needs("--start_bp")->needs("--end_bp");
    app_.get_formatter()->column_width(80);

}

void
DesignRNAScaffold::parse_command_line (
    int argc,
    const char **argv) {
}

void
DesignRNAScaffold::run () {
    // sets up all variables required
    setup();

    LOG_INFO << "###################";
    LOG_INFO << "# starting search #";
    LOG_INFO << "###################";
    sol_info_ = SolutionInfo();
    sol_info_.design_num = 0;
    // main run loop
    while (true) {
        LOG_DEBUG << "*** starting motif search ***";
        auto mg_w_sol = _get_motif_graph_solution();
        if (mg_w_sol == nullptr) {
            LOG_DEBUG << "no more solutions generated by motif search";
            break;
        }
        if (parameters_.seq_opt.skip || parameters_.search.only_tether_opt) {
            _record_solution(*mg_w_sol);
            LOG_DEBUG << "*** finished design: " << sol_info_.design_num << " ***";
            sol_info_.design_num += 1;
            if (parameters_.core.designs <= sol_info_.design_num) {
                LOG_DEBUG << " successfully finished main design loop!";
                break;
            }
            continue;
        }
        // remap node indexes with new base pair step motifs
        auto bp_step_indexes = GraphIndexes();
        _get_graph_indexes_after_bp_steps(*mg_w_sol, starting_indexes_, bp_step_indexes);
        LOG_DEBUG << "*** starting sequence optimization ***";
        sol_info_.seqeunce_opt_num = 0;
        auto attempts = 1;
        auto seq_opt_fail = false;
        for (int i = 0; i < parameters_.seq_opt.sequences_per_design; i++) {
            LOG_DEBUG << "starting sequence optimization attempt: " << i;
            auto mg_seq_opt = _perform_sequence_opt(*mg_w_sol, bp_step_indexes);
            if (mg_seq_opt == nullptr) {
                LOG_DEBUG << "did not find a valid solution for sequence opt";
                continue;
            }
            if (sol_info_.sequence_opt_score > 7) {
                attempts += 1;
                if (attempts > 3) {
                    LOG_DEBUG << "no viable sequence determined skipping this design ";
                    seq_opt_fail = true;
                    break;
                }
                LOG_DEBUG << "seq opt did not a good enough solution: " << sol_info_.sequence_opt_score << " redoing attempt " << attempts;
                i -= 1;
                continue;
            }
            sol_info_.seqeunce_opt_num += 1;
            if (!parameters_.thermo_fluc.perform) {
                _record_solution(*mg_seq_opt);
                continue;
            }
            auto mg_thermo = _perform_thermo_fluc_sim(*mg_seq_opt, bp_step_indexes);
            _record_solution(*mg_thermo);

        }
        if (seq_opt_fail) {
            continue;
        }
        LOG_DEBUG << "*** finished design: " << sol_info_.design_num << " ***";
        sol_info_.design_num += 1;
        if (parameters_.core.designs <= sol_info_.design_num) {
            LOG_DEBUG << " successfully finished main design loop!";
            break;
        }

    }

    if (parameters_.core.designs <= sol_info_.design_num) {
        LOG_INFO << "found " << sol_info_.design_num
                 << " designs; if you want more please use --designs num_of_design";

    }
    else {
        LOG_INFO << "found " << sol_info_.design_num << " designs; did not reach desired num of "
                 << parameters_.core.designs;

    }
    LOG_DEBUG << "NON-ERROR EXIT";
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// setup functions
////////////////////////////////////////////////////////////////////////////////////////////////////

void
DesignRNAScaffold::setup () {
    LOG_INFO << "#########";
    LOG_INFO << "# setup #";
    LOG_INFO << "#########";

    if (!parameters_.io.new_ensembles_file.empty()) {
        _build_new_ensembles(parameters_.io.new_ensembles_file);
    }
    // setups up start_, end_ and msg_;
    _setup_from_pdb();
    // extra pdbs that might be included in build path
    _setup_extra_pdbs();

    mg_ = msg_->to_motif_graph();
    mg_->set_option_value("sterics", false);
    mg_->increase_level();

    // sets up outputing files
    LOG_INFO << "out file that contains all build solutions -> " + parameters_.io.out_file
            + " can be set with --out_file";
    LOG_INFO << "score file that contains summary information about solutions -> "
            + parameters_.io.score_file + " can be set with --score_file";
    out_.open(parameters_.io.out_file);
    score_out_.open(parameters_.io.score_file);

    // keeps track of what is being recorded
    score_out_ << SolutionInfo::col_names() << std::endl;

    // setup search
    search_ = _setup_search();
    problem_ = _setup_problem();
    search_->setup(problem_);

    //sequence optimziation setup
    seq_optimizer_ = std::make_shared<sequence_optimization::SequenceOptimizer3D>();
    seq_optimizer_->set_option_value("steps", parameters_.seq_opt.steps);

    if (!parameters_.seq_opt.skip) {
        LOG_INFO
                << "sequence optimization cutoff -> " << seq_optimizer_->get_float_option("cutoff");
        LOG_INFO << "sequence optimization steps -> " << seq_optimizer_->get_int_option("steps");
    }
    else {
        LOG_INFO << "set --skip_sequence_optimization -> skipping sequence optimization step!";
        LOG_WARNING << "do not trust the sequence in the outputed PDBs!";
    }

    //thermo sim setup
    auto thermo_scorer = std::make_shared<thermo_fluctuation::graph::OldFrameScorer>();
    auto sterics = std::make_shared<thermo_fluctuation::graph::sterics::NoSterics>();
    thermo_sim_ = std::make_shared<thermo_fluctuation::graph::Simulation>(thermo_scorer, sterics);

}

void
DesignRNAScaffold::_setup_from_pdb () {
    auto struc = rm_.get_structure(parameters_.core.pdb, "scaffold");
    LOG_INFO << "loaded pdb from file: " << parameters_.core.pdb;
    if (struc->ends().empty()) {
        LOG_ERROR << "pdb has no ends available to build off of";
        exit(0);
    }
    auto end_names = String();
    for (auto const &end : struc->ends()) {
        end_names += end->name() + " ";
    }
    LOG_INFO << "available ends to build off: " << end_names;
    _check_bp(parameters_.core.start_bp, struc, "start");
    _check_bp(parameters_.core.end_bp, struc, "end");

    // TODO allow for construction of motif from RNA structure to allow for changes in base pair types
    rm_.add_motif(parameters_.core.pdb, "scaffold", util::MotifType::TCONTACT);
    auto m = rm_.motif("scaffold", "", parameters_.core.end_bp);
    int ei1, ei2;

    try {
        ei1 = m->get_end_index(parameters_.core.start_bp);
    }
    catch (structure::RNAStructureException const &e) {
        LOG_ERROR << "starting basepair: " << parameters_.core.start_bp << " was NOT FOUND!";
        exit(0);
    }
    LOG_INFO << "start basepair: " << parameters_.core.start_bp << " is found!";
    try {
        ei2 = m->get_end_index(parameters_.core.end_bp);
    }
    catch (structure::RNAStructureException const &e) {
        LOG_ERROR << "end basepair: " << parameters_.core.end_bp << " was NOT FOUND!";
        exit(0);
    }
    LOG_INFO << "end basepair: " << parameters_.core.end_bp << " is found!";

    // TODO not a great way of doing this ...
    lookup_ = std::make_shared<util::StericLookupNew>();
    auto total_beads = 0;
    auto protein_beads = 0;
    if (!parameters_.search.no_sterics) {
        protein_beads += m->protein_beads().size();
        for(auto const & b : m->beads()) {
            if (b.btype() == structure::BeadType::PHOS) { continue; }
            lookup_->add_point(b.center());
            total_beads += 1;
        }
        for(auto const & b : m->protein_beads()) {
            lookup_->add_point(b.center());
        }
        LOG_INFO << total_beads
                 << " rna steric beads were found to block residue overlap with supplied structure";
        LOG_INFO << protein_beads
                 << " protein steric beads were found to block residue overlap with supplied structure";
    }
    else {
        LOG_INFO << "no steric beads will be added based on --no_sterics flag";
    }

    if(parameters_.search.only_tether_opt) {
        LOG_INFO << "outputing scaffold.pdb for reference to align designs";
        m->to_pdb("scaffold.pdb", 1, 1);
    }

    starting_indexes_.start = data_structure::NodeIndexandEdge{0, ei1};
    starting_indexes_.end = data_structure::NodeIndexandEdge{0, ei2};

    msg_ = std::make_shared<motif_data_structure::MotifStateGraph>();
    msg_->add_state(m->get_state());
    if (!parameters_.search.starting_helix.empty()) {
        auto h = rm_.motif_state(parameters_.search.starting_helix);
        msg_->add_state(h, starting_indexes_.start.node_index, starting_indexes_.start.edge_index);
        starting_indexes_.start = data_structure::NodeIndexandEdge{1, 1};
    }
}

void
DesignRNAScaffold::_setup_extra_pdbs() {
    auto pdb_files = base::split_str_by_delimiter(parameters_.core.extra_pdbs, ",");
    if(!pdb_files.empty()) {
        LOG_INFO << " loading pdbs from --extra_pdb flag: " << parameters_.core.extra_pdbs;
    }
    for(auto const & pdb_file : pdb_files) {
        if(!base::file_exists(pdb_file)) {
            LOG_ERROR << "invalid pdb path: " << pdb_file;
            exit(0);
        }
        LOG_INFO << "loading " << pdb_file;
        rm_.add_motif(pdb_file);
    }
}

void
DesignRNAScaffold::_check_bp (
    String const &name,
    structure::RNAStructureOP const &struc,
    String const &type) const {

    auto bps = structure::BasepairOPs();
    try {
        bps = struc->get_basepair(name);
    }
    catch (std::runtime_error const &error) {
        LOG_ERROR << error.what();
        LOG_ERROR << "cannot find " + type + " basepair " + name;
        exit(0);
    }

    if (bps[0]->bp_type() != "cW-W" && parameters_.search.no_basepair_checks) {
        bps[0]->bp_type("cW-W");
        LOG_WARNING << "basepair " + name
                + " is not a watson and crick bp, but forcing it to be, this may be bad";
    }

    if (bps[0]->bp_type() != "cW-W") {
        LOG_ERROR << "basepair " + name + " is not a watson and crick bp ";
        exit(0);
    }

}

motif_search::SearchOP
DesignRNAScaffold::_setup_search () {
    auto search = motif_search::SearchOP(nullptr);
    auto filter = _setup_sol_filter(parameters_.search.solution_filter);
    LOG_INFO << "search type -> " + parameters_.search.type;
    LOG_INFO << "search solution filter -> " + parameters_.search.solution_filter;
    if (parameters_.search.solution_filter == "NoFilter") {
        LOG_WARNING << "no solution filter is set, this will lead to duplicate solutions under "
                    << "certain conditions";
    }

    if (parameters_.search.type == "path_finding") {
        using namespace motif_search::path_finding;
        auto scorer = std::make_shared<AstarScorer>();
        auto selector = default_selector();
        auto pf_search = std::make_shared<Search>(scorer, selector, filter);
        pf_search->set_option_value("max_node_level", parameters_.search.max_motifs);
        search = motif_search::SearchOP(pf_search->clone());
    }

    if (parameters_.search.type == "exhaustive") {
        using namespace motif_search::exhaustive;
        auto score_factory = ScorerFactory();
        auto scorer = score_factory.get_scorer(parameters_.search.exhaustive_scorer);
        //TODO setup default sol topology if user does not specify motif_path
        auto sol_template = std::make_shared<motif_search::SolutionTopologyTemplate>();
        if (!parameters_.search.motif_path.empty() /* added by CJ 08/20 parameters_.motif_path != ""*/) {
            sol_template = _setup_sol_template_from_path(parameters_.search.motif_path);
        }
        auto factory = motif_search::SolutionToplogyFactory();
        auto sol_toplogy = factory.generate_toplogy(*sol_template);
        auto e_search = std::make_shared<Search>(scorer, *sol_toplogy, filter);
        search = motif_search::SearchOP(e_search->clone());
    }

    if (parameters_.search.type == "mc") {
        using namespace motif_search::monte_carlo;
        auto scorer = ScorerOP(nullptr);
        if (parameters_.search.mc_scorer == "default") {
            scorer = std::make_shared<DefaultScorer>();
        }
        else if (parameters_.search.mc_scorer == "scaled_scorer") {
            scorer = std::make_shared<ScaledScorer>(
                parameters_.search.scaled_score_d,
                parameters_.search.scaled_score_r);
        }
        if (parameters_.search.motif_path.empty()) {
            LOG_ERROR << "must include a motif path when using Monte Carlo search";
            exit(0);
        }
        auto sol_template = _setup_sol_template_from_path(parameters_.search.motif_path);
        auto factory = motif_search::SolutionToplogyFactory();
        factory.set_option_value("max_helix_size", parameters_.search.max_helix_length);
        factory.set_option_value("min_helix_size", parameters_.search.min_helix_length);
        auto sol_toplogy = factory.generate_toplogy(*sol_template);
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
DesignRNAScaffold::_setup_problem () {
    auto start = starting_indexes_.start;
    auto end = starting_indexes_.end;
    auto start_bp = msg_->get_node(start.node_index)->data()->get_end_state(start.edge_index);
    auto end_bp = msg_->get_node(end.node_index)->data()->get_end_state(end.edge_index);

    bool target_an_aligned_end = false;
    if (end.edge_index == msg_->get_node(end.node_index)->data()->block_end_add()) {
        target_an_aligned_end = true;
    }
    return std::make_shared<motif_search::Problem>(start_bp, end_bp, lookup_, target_an_aligned_end);

}

motif_search::SolutionTopologyTemplateOP
DesignRNAScaffold::_setup_sol_template_from_path (
    String const &motif_path) {
    auto sol_template = std::make_shared<motif_search::SolutionTopologyTemplate>();
    auto spl = base::split_str_by_delimiter(motif_path, ",");
    int i = 0;
    auto lib_names = resources::MotifStateSqliteLibrary::get_libnames();
    LOG_DEBUG << "building motif path from string: " << motif_path;
    for (auto const &e : spl) {
        if (lib_names.find(e) != lib_names.end()) {  // is a library
            LOG_DEBUG << e << " is determined to be a motif library";
            if (i == 0) {
                sol_template->add_library(e);
            }
            else {
                sol_template->add_library(e, data_structure::NodeIndexandEdge{i - 1, 1});
            }
        }

        else if (e == "tc_hairpin_hairpin") {
            LOG_ERROR << "tc_hairpin_hairpin not supported";
            exit(0);
        }
        // found user specified ensemble
        else if (new_motif_ensembles_.find(e) != new_motif_ensembles_.end()) {
            LOG_DEBUG << e << " is determined to be a motif library";
            if (i == 0) {
                sol_template->add_ensemble(new_motif_ensembles_[e]);
            }
            else {
                sol_template->add_ensemble(new_motif_ensembles_[e],
                                           data_structure::NodeIndexandEdge{i - 1, 1});
            }
        }
        else {
            // TODO should try and get motif from every end!
            auto ms = motif::MotifStateOP(nullptr);
            try {
                ms = rm_.motif_state(e);
                // Need to do this in case I am using two copies of the same motif
                ms->new_uuids();
                //std::cout << ms->uuid() << std::endl;
            }
            catch(resources::ResourceManagerException const & error) {
                LOG_ERROR << "unclear what " << e << " is it is not recognized as a motif library," <<
                        " ensemble or supplied motif";
                exit(0);
            }
            LOG_DEBUG << e << " is determined to be a supplied motif";
            if (i == 0) {
                sol_template->add_motif_state(ms);
            }
            else {
                sol_template->add_motif_state(ms, data_structure::NodeIndexandEdge{i - 1, 1});
            }
        }

        i++;
    }

    return sol_template;
}

motif_search::SolutionFilterOP
DesignRNAScaffold::_setup_sol_filter (
    String const &solution_filter_name) {
    auto sf = motif_search::SolutionFilterOP(nullptr);
    if (solution_filter_name == "NoFilter") {
        auto new_sf = std::make_shared<motif_search::NoExclusionFilter>();
        sf = motif_search::SolutionFilterOP(new_sf->clone());
    }
    else if (solution_filter_name == "RemoveDuplicateHelices") {
        auto new_sf = std::make_shared<motif_search::RemoveDuplicateHelices>();
        sf = motif_search::SolutionFilterOP(new_sf->clone());
    }
    else {
        LOG_ERROR << "not a valid solution filter name: " << solution_filter_name;
        exit(0);
    }

    return sf;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// run functions
////////////////////////////////////////////////////////////////////////////////////////////////////

motif_data_structure::MotifGraphOP
DesignRNAScaffold::_get_motif_graph_solution() {
    auto solution = search_->next();
    // search could not find a solution
    if (solution == nullptr) {
        return nullptr;
    }
    auto sol_mg = solution->graph->to_motif_graph();
    // for some reason solution is size 0 this reall shouldnt happen
    if (sol_mg->size() == 0) {
        return nullptr;
    }
    // load info into solution_info
    sol_info_.design_score = solution->score;
    _get_motif_names(sol_mg);

    if (sol_info_.design_num < 10) {
        LOG_INFO << "found a solution: " << sol_info_.motif_names << " with score: "
                 << sol_info_.design_score;
    }
    else if (sol_info_.design_num % 10 == 0) {
        LOG_INFO << "found " << sol_info_.design_num << " solutions ";
    }

    //msg_->to_motif_graph()->to_pdb("scaffold.pdb", 1, 1);
    //sol_mg->to_pdb("tether.pdb", 1, 1x);
    //exit(0);

    auto mg_w_sol = motif_data_structure::MotifGraphOP(nullptr);
    if(parameters_.search.only_tether_opt) {
        mg_w_sol = std::make_shared<motif_data_structure::MotifGraph>(*sol_mg);
    }
    else {
        auto start = starting_indexes_.start;
        auto end = starting_indexes_.end;
        mg_->add_motif_graph(*sol_mg, start.node_index, start.edge_index);
        mg_w_sol = std::make_shared<motif_data_structure::MotifGraph>(*mg_);
        mg_w_sol->add_connection(mg_->last_node()->index(), end.node_index,
                                 mg_->last_node()->data()->end_name(1),
                                 mg_->get_node(end.node_index)->data()->end_name(end.edge_index));
        _fix_flex_helices_mtype(mg_w_sol);

        mg_->remove_level(1);
    }
    // must output the pdb before replacing idealized helices
    if (parameters_.io.dump_pdbs && parameters_.seq_opt.skip) {
        mg_w_sol->to_pdb("design." + std::to_string(sol_info_.design_num) + ".pdb", 1, 1);
    }
    // currently necessary to generate correct designable secondary structure
    mg_w_sol->replace_ideal_helices();
    LOG_DEBUG << "potential solution: #" << sol_info_.design_num;
    auto ss = mg_w_sol->designable_secondary_structure();
    sol_info_.designable_sequence = ss->sequence();
    sol_info_.dot_bracket = ss->dot_bracket();
    LOG_DEBUG << ss->sequence();
    LOG_DEBUG << ss->dot_bracket();
    return mg_w_sol;
}

void
DesignRNAScaffold::_get_graph_indexes_after_bp_steps (
    motif_data_structure::MotifGraph const &mg,
    GraphIndexes const &starting_indexes,
    GraphIndexes &new_indexes /* return */) {
    auto end = starting_indexes.end;
    auto last_m =
        mg.get_node(end.node_index)->connections()[end.edge_index]->partner(end.node_index);
    new_indexes.start = data_structure::NodeIndexandEdge{last_m->index(), 1};
    new_indexes.end = end;
}

motif_data_structure::MotifGraphOP
DesignRNAScaffold::_perform_sequence_opt (
    motif_data_structure::MotifGraph const &mg_w_sol,
    GraphIndexes const &bp_step_indexes) {
    auto opt_seq_scorer = std::make_shared<sequence_optimization::InternalTargetScorer>(
        bp_step_indexes.start.node_index, bp_step_indexes.start.edge_index,
        bp_step_indexes.end.node_index, bp_step_indexes.end.edge_index,
        problem_->target_an_aligned_end);
    auto mg_seq_opt = std::make_shared<motif_data_structure::MotifGraph>(mg_w_sol);
    auto sols = seq_optimizer_->get_optimized_sequences(mg_seq_opt, opt_seq_scorer);
    if (sols.empty()) {
        LOG_DEBUG << "No viable solutions sequence solutions founds!";
        return nullptr;
    }
    mg_seq_opt->replace_helical_sequence(sols[0]->sequence);
    sol_info_.sequence = sols[0]->sequence;
    sol_info_.sequence_opt_score = sols[0]->dist_score;
    return mg_seq_opt;
}

motif_data_structure::MotifGraphOP
DesignRNAScaffold::_perform_thermo_fluc_sim (
    motif_data_structure::MotifGraph &mg,
    GraphIndexes const &bp_step_indexes) {
    LOG_DEBUG << "*** starting thermo fluc sim ***";
    auto index_hash = std::map<int, int>();
    auto mseg = motif_data_structure::MotifStateEnsembleGraph();
    _get_mseg(mg, mseg, index_hash);
    auto tf_indexes = GraphIndexes();
    tf_indexes.start =
        data_structure::NodeIndexandEdge{index_hash[bp_step_indexes.start.node_index],
                                         bp_step_indexes.start.edge_index};
    tf_indexes.end = data_structure::NodeIndexandEdge{index_hash[bp_step_indexes.end.node_index],
                                                      bp_step_indexes.end.edge_index};
    thermo_sim_->setup(mseg, tf_indexes.start, tf_indexes.end);
    thermo_sim_->next();
    auto count = 0;
    auto best = thermo_sim_->get_score();
    auto best_mg = thermo_sim_->get_motif_graph();
    auto under_cutoff = 0;
    for (int s = 0; s < parameters_.thermo_fluc.steps; s++) {
        under_cutoff = thermo_sim_->next();
        if (under_cutoff) {
            count += 1;
        }
        if (thermo_sim_->get_score() < best) {
            LOG_VERBOSE << "best dist score: " << best;
            best = thermo_sim_->get_score();
            best_mg = thermo_sim_->get_motif_graph();
        }
    }
    // connection is not perserved through simulation ... bring it back
    best_mg->add_connection(
        tf_indexes.start.node_index,
        tf_indexes.end.node_index,
        best_mg->get_node(tf_indexes.start.node_index)->data()->ends()[tf_indexes.start.edge_index]->name(),
        best_mg->get_node(tf_indexes.end.node_index)->data()->ends()[tf_indexes.end.edge_index]->name());

    sol_info_.thermo_fluc_best_score = best;
    sol_info_.thermo_fluc_hits = count;
    LOG_DEBUG << "best dist score: " << best;
    LOG_DEBUG << "hits: " << count << ", %: = " << ((float) count / (float) 1000000) * 100.0f;
    LOG_DEBUG << "*** finished thermo fluc sim ***";
    return best_mg;
}

void
DesignRNAScaffold::_record_solution (
    motif_data_structure::MotifGraph &mg) {

    if (parameters_.io.dump_pdbs && !parameters_.seq_opt.skip) {
        mg.to_pdb("design." + std::to_string(sol_info_.design_num) + ".pdb", 1, 1);
    }
    score_out_ << sol_info_ << std::endl;
    out_ << mg.to_str() << std::endl;
}

void
DesignRNAScaffold::_fix_flex_helices_mtype (
    motif_data_structure::MotifGraphOP mg) {
    // Some versions of RNAMake have this issue, that flex helices arent set as "helices"
    // this can raise merger issues
    for (auto &n : *mg) {
        if (n->data()->name().substr(0, 5) == "HELIX") {
            n->data()->mtype(util::MotifType::HELIX);
        }
    }
}

void
DesignRNAScaffold::_get_motif_names (
    motif_data_structure::MotifGraphOP mg) {
    sol_info_.motif_names = "";
    for (auto const &n : *mg) {
        sol_info_.motif_names += n->data()->name() + ";";
    }
}

void
DesignRNAScaffold::_build_new_ensembles (
    String const &new_ensemble_path) {
    auto lines = base::get_lines_from_file(new_ensemble_path);
    auto ensemble = motif::MotifStateEnsembleOP(nullptr);
    auto ms = motif::MotifStateOP(nullptr);
    auto all_ms = motif::MotifStateOPs();
    auto scores = Floats();
    for (auto const &l : lines) {
        auto spl = base::split_str_by_delimiter(l, " ");
        if (spl.size() != 3) {
            continue;
        }
        all_ms.resize(0);
        scores.resize(0);
        for (auto const &m_name : base::split_str_by_delimiter(spl[2], ",")) {
            ms = rm_.motif_state(m_name, "", "");
            for (auto const &end_name : ms->end_names()) {
                try {
                    all_ms.push_back(rm_.motif_state(m_name, "", end_name));
                    scores.push_back(0);
                }
                catch (std::runtime_error const &error) {
                    LOGW << "Error: " << error.what();
                    continue;
                }
            }
        }
        new_motif_ensembles_[spl[0]] = std::make_shared<motif::MotifStateEnsemble>(all_ms, scores);
    }

}

void
DesignRNAScaffold::_get_mseg (
    motif_data_structure::MotifGraph &mg,
    motif_data_structure::MotifStateEnsembleGraph &mseg /* return */,
    std::map<int, int> &index_hash) {
    int i = 0, j = 0;
    for (auto const &n : mg) {
        // build / get motif state ensemble
        auto mse = motif::MotifStateEnsembleOP(nullptr);

        if (n->data()->mtype() == util::MotifType::HELIX) {
            // is not a basepair step
            if (n->data()->residues().size() > 4) {
                LOG_ERROR << "supplied a helix motif: " + n->data()->name()
                        + " that is not a basepair "
                          << "step, this is no supported";
                exit(0);
            }

            try {
                mse = rm_.motif_state_ensemble(n->data()->end_ids()[0]);
            }
            catch (resources::ResourceManagerException const &e) {
                LOG_ERROR << "cannot find motif state ensemble for basepair with id: "
                        + n->data()->end_ids()[0] <<
                          "check to make sure its a Watson-Crick basepair step";
                exit(0);
            }
        }
        else {
            mse = std::make_shared<motif::MotifStateEnsemble>(n->data()->get_state());
        }

        if (i == 0) {
            j = mseg.add_ensemble(*mse);
        }
        else {
            int pi = index_hash[mg.parent_index(n->index())];
            int pie = mg.parent_end_index(n->index());
            j = mseg.add_ensemble(*mse, data_structure::NodeIndexandEdge{pi, pie});
        }

        index_hash[n->index()] = j;
        i++;

    }

}

// There used to be a main function here, is its absence correct?
