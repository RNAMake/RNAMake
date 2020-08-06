//
// Created by Joseph Yesselman on 2/23/19.
//

#include <algorithm>

#include "base/backtrace.h"
#include "base/log.h"
#include "resources/resource_manager.h"
#include "motif_data_structure/motif_state_graph.hpp"
#include "motif_search/motif_state_monte_carlo.h"
#include "sequence_optimization/sequence_optimizer_3d.hpp"
#include "apt_new_interface/apt_new_interface.h"

AptNewInterface::AptNewInterface():
    rm_(resources::Manager::instance()),
    rng_(util::RandomNumberGenerator()) {}

void
AptNewInterface::setup_options() {
    add_option("scaffold", String(""), base::OptionType::STRING, true);
    add_option("docked_motif", String(""), base::OptionType::STRING, true);
    //add_option("scaffold_end", String(""), base::OptionType::STRING, true);
    add_option("connections", String(""), base::OptionType::STRING, true);

    // general options
    add_option("out_file", "default.out", base::OptionType::STRING, false);
    add_option("score_file", "default.scores", base::OptionType::STRING, false);
    add_option("designs", 1, base::OptionType::INT, false);

    add_option("search_cutoff_1", 5.0f, base::OptionType::FLOAT, false);
    add_option("search_cutoff_2", 5.0f, base::OptionType::FLOAT, false);

    add_option("motifs_1", 2, base::OptionType::INT, false);
    add_option("motifs_2", 3, base::OptionType::INT, false);

    add_option("rounds", 50, base::OptionType::INT, false);

    add_option("testing", false, base::OptionType::BOOL, false);


}

void
AptNewInterface::parse_command_line(
        int argc,
        const char **argv) {

    base::Application::parse_command_line(argc, argv);
}

void
AptNewInterface::run() {
    _setup_new_motifs();

    std::ofstream out, out_score;
    out.open("default.out");
    out_score.open("default.scores");
    out_score << "design_num,score_1,score_2,size,motifs" << std::endl;

    auto ms_libraries_1 = _get_libraries(get_int_option("motifs_1"));
    auto ms_libraries_2 = _get_libraries(get_int_option("motifs_2"));
    auto problem = _get_design_problem();
    auto org_motif_num = problem->msg->size();
    auto org_res = 0;
    for(auto const & n : *problem->msg) {
        org_res += n->data()->cur_state->size();
    }
    std::cout << org_res << std::endl;
    exit(0);
    auto design_num = 0;

    for(int i = 0; i < 1; i++) {
        lookup_ = util::StericLookupNew();
        _setup_sterics(problem->msg);

        auto search_1 = _setup_search(problem->path_1, lookup_, problem->msg, ms_libraries_1,
                                      get_float_option("search_cutoff_1"));

        auto sol = search_1->next_state();
        if(sol == nullptr) { continue; }
        lookup_ = util::StericLookupNew();
        auto mg = sol->msg->to_motif_graph();
        for (int j = 0; j < mg->size(); j++) {
            for(auto const & b : mg->get_node(j)->data()->beads()) { lookup_.add_point(b.center()); }
        }

        for(int j = 0; j < 10; j++) {
            auto search_2 = _setup_search(problem->path_2, lookup_, sol->msg, ms_libraries_2,
                                          get_float_option("search_cutoff_2"));
            auto sol_2 = search_2->next();

            if(sol_2 == nullptr) { continue; }
            auto res_count = _get_residue_count(sol_2->mg) - org_res;

            out << sol_2->mg->to_str() << std::endl;
            out_score << design_num << "," << sol->score << "," << sol_2->score << "," << res_count << ",";
            out_score << _get_motif_names(sol_2->mg) << std::endl;
            design_num += 1;

        }
    }

    out.close();
    out_score.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<motif::MotifStateOPs>
AptNewInterface::_get_libraries(
        int motif_num) {

    auto num = motif_num*2 + 1;
    auto motif_lib_names = Strings(num);
    auto motif_types = std::vector<util::MotifType>(num);
    for(int i = 0; i < num; i++) {
        if (i % 2 == 0) {
            //motif_lib_names[i] = parameters_.helices;
            motif_lib_names[i] = "flex_helices";
            motif_types[i]     = util::MotifType::HELIX;
        }
        else            {
            motif_lib_names[i] = "twoway";
            motif_types[i]     = util::MotifType::TWOWAY;
        }
    }

    auto libraries = std::vector<motif::MotifStateOPs>();
    auto motif_states = motif::MotifStateOPs();
    int i = 0;
    for(auto const & name : motif_lib_names) {
        auto ms_lib =  resources::MotifStateSqliteLibrary(name);
        ms_lib.load_all();
        motif_states = motif::MotifStateOPs();

        if(motif_types[i] == util::MotifType::HELIX) {
            for (auto const & ms : ms_lib) {
                //if(ms->size() >= parameters_.min_helix_size && ms->size() <= parameters_.max_helix_size) {
                motif_states.push_back(ms);
                //}
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

void
AptNewInterface::_setup_sterics(
        motif_data_structure::MotifStateGraphOP msg) {
    auto beads = math::Points();
    for (auto & n : *msg) {
        for(auto const & b : n->data()->cur_state->beads()) {
            beads.push_back(b);
        }
    }
    lookup_.add_points(beads);
}


void
AptNewInterface::_setup_new_motifs() {
    auto mf = motif::MotifFactory();

    // add new motifs to resource manager
    auto scaffold_rm = rm_.get_structure(get_string_option("scaffold"), "scaffold", 3);
    auto scaffold_m = std::make_shared<motif::Motif>(*scaffold_rm);
    scaffold_m->name("scaffold");
    scaffold_m->mtype(util::MotifType::TCONTACT);
    scaffold_m->score(0);
    //scaffold_m->block_end_add(0);
    mf._setup_secondary_structure(scaffold_m);
    rm_.register_motif(scaffold_m);

    auto rs = rm_.get_structure(get_string_option("docked_motif"), "docked_motif");

    auto docked_motif = std::make_shared<motif::Motif>(*rs);
    docked_motif->mtype(util::MotifType::HAIRPIN);
    docked_motif->block_end_add(0);
    docked_motif->score(0);
    mf._setup_secondary_structure(docked_motif);

    rm_.register_motif(docked_motif);
}

/*AptNewInterfaceProblemOP
AptNewInterface::_get_design_problem() {
    auto prna = rm_.motif("prna", "", "");
    auto end_names = Strings();
    for(auto const & end : prna->ends()) {
        end_names.push_back(end->name());
        std::cout << end->name() << " ";
    }
    std::cout << std::endl;
    exit(0);

    std::random_device rd;
    std::mt19937 g(rd());

    std::shuffle(end_names.begin(), end_names.end(), g);
    prna = rm_.motif("prna", "", end_names[0]);

    auto scaffold = rm_.motif("scaffold", "", "A1-A129");
    auto docked_motif = rm_.motif("docked_motif", "", "");

    auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
    msg->set_option_value("sterics", false);
    msg->add_state(scaffold->get_state());
    msg->add_state(rm_.motif_state("HELIX.IDEAL.1"), 0, 2);
    msg->add_state(prna->get_state());
    msg->add_state(docked_motif->get_state(), -1, -1, 1);
    msg->increase_level();

    auto node_infos = std::vector<NodeIndexandEdge>(4);
    node_infos[0] = NodeIndexandEdge{2, prna->get_end_index(end_names[1])};
    node_infos[1] = NodeIndexandEdge{2, prna->get_end_index(end_names[2])};
    node_infos[2] = NodeIndexandEdge{0, scaffold->get_end_index("A41-A87")};
    node_infos[3] = NodeIndexandEdge{3, 0}; // docked hairpin

    auto pairs = std::vector<Ints> {
            Ints{0,2}, Ints{2,0}, Ints{0,3}, Ints{1,2}, Ints{2,1}, Ints{1,3}};

    std::shuffle(pairs.begin(), pairs.end(), g);

    auto connections = std::vector<ConnectionInfo>();
    auto used = std::map<int, int>();

    for(auto const & p : pairs) {
        if(used.find(p[0]) != used.end() || used.find(p[1]) != used.end()) { continue; }
        connections.push_back(ConnectionInfo{node_infos[p[0]], node_infos[p[1]]});
        if(connections.size() == 2) { break; }
    }

    auto p = std::make_shared<AptNewInterfaceProblem>(msg, connections[0], connections[1]);
    return p;

}*/

AptNewInterfaceProblemOP
AptNewInterface::_get_design_problem() {
    auto connection_str = get_string_option("connections");
    auto available_bps = StringIntMap{{"A7-C10", 2}, {"B14-C7", 2}, {"A10-B9", 2}, {"A60-A65", 3}, {"A41-A87", 0}};
    auto seen_bps = StringIntMap();
    auto scaffold_build_bp_name = String("A59-A68");

    auto connections = std::vector<std::vector<std::pair<String, int>>>();

    auto spl = base::split_str_by_delimiter(connection_str, ";");
    if(spl.size() != 2) { LOG_ERROR << "must supply two connections seperated by a ;"; }

    for(auto const & s : spl) {
        auto bp_spl = base::split_str_by_delimiter(s, ":");
        auto connection = std::vector<std::pair<String, int>>();
        for(auto const bp_str : bp_spl) {
            if(available_bps.find(bp_str) == available_bps.end()) { LOG_ERROR << bp_str << "is not a valid end"; }
            seen_bps[bp_str] = 1;
            auto pair = std::pair<String, int>(bp_str, available_bps[bp_str]);
            connection.push_back(pair);
        }
        connections.push_back(connection);
    }

    auto prna_ends = Strings{"A7-C10", "B14-C7", "A10-B9"};
    auto start_prna_end = String();
    for(auto const & end : prna_ends) {
        if(seen_bps.find(end) != seen_bps.end()) { continue; }
        start_prna_end = end;
    }

    auto scaffold = rm_.motif("scaffold", "", "A1-A129");
    auto prna = rm_.motif("prna", "", start_prna_end);
    auto docked_motif = rm_.motif("docked_motif", "", "");

    auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
    msg->set_option_value("sterics", false);
    msg->add_state(scaffold->get_state());
    msg->add_state(rm_.motif_state("HELIX.IDEAL.3"), 0, 2);
    msg->add_state(prna->get_state());
    msg->add_state(docked_motif->get_state(), -1, -1, 1);
    msg->increase_level();

    auto connection_infos = std::vector<ConnectionInfo>();
    for(auto const & c : connections) {
        auto start = NodeIndexandEdge{c[0].second, msg->get_node(c[0].second)->data()->get_end_index(c[0].first)};
        auto end   = NodeIndexandEdge{c[1].second, msg->get_node(c[1].second)->data()->get_end_index(c[1].first)};
        connection_infos.push_back(ConnectionInfo{start, end});
    }

    auto p = std::make_shared<AptNewInterfaceProblem>(msg, connection_infos[0], connection_infos[1]);
    return p;
}

motif_data_structure::MotifStateGraphOP
AptNewInterface::_setup_graph() {

    //auto prna = rm_.motif("prna", "", "A7-C10");

    auto scaffold = rm_.motif("scaffold", "", "A1-A129");
    auto prna = rm_.motif("prna", "", "");
    auto docked_motif = rm_.motif("docked_motif", "", "");

    auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
    msg->set_option_value("sterics", false);
    msg->add_state(scaffold->get_state());
    msg->add_state(rm_.motif_state("HELIX.IDEAL.1"), 0, 2);
    msg->add_state(prna->get_state());
    msg->add_state(docked_motif->get_state(), -1, -1, 1);
    msg->increase_level();

    return msg;

}

motif_search::MotifStateMonteCarloOP
AptNewInterface::_setup_search(
        ConnectionInfo const & ci,
        util::StericLookupNew const & sl,
        motif_data_structure::MotifStateGraphOP msg,
        std::vector<motif::MotifStateOPs> const & ms_libraries,
        float cutoff) {

    bool target_an_aligned_end = false;
    if(ci.end.ei == msg->get_node(ci.end.ni)->data()->block_end_add()) {
        target_an_aligned_end = true;
    }

    auto mc = std::make_shared<motif_search::MotifStateMonteCarlo>(ms_libraries);
    mc->setup(msg, ci.start.ni, ci.end.ni, ci.start.ei, ci.end.ei, target_an_aligned_end);
    mc->lookup(sl);
    mc->set_option_value("accept_score", cutoff);

    return mc;
}

int
AptNewInterface::_find_min_motifs_per_path(
        ConnectionInfo const & ci,
        util::StericLookupNew const & sl,
        motif_data_structure::MotifStateGraphOP msg,
        float cutoff,
        int motif_num) {

    auto done = false;
    while(!done) {
        auto ms_libraries = _get_libraries(motif_num);
        auto search = _setup_search(ci, sl, msg, ms_libraries, cutoff);
        if(get_bool_option("testing")) { search->set_option_value("stages", 10); }
        auto sol = search->next();

        if(sol == nullptr) { return motif_num+1; }
        LOGI << "solution found with " << motif_num << " motifs";
        if(motif_num == 1) { return motif_num; }
        motif_num -= 1;
    }
    return motif_num;

}


int
AptNewInterface::_get_residue_count(
        motif_data_structure::MotifGraphOP mg) {
    auto res_count = 0;
    for(auto const & n : *mg) {
        res_count += n->data()->residues().size();
    }
    return res_count;
}

String
AptNewInterface::_get_motif_names(
        motif_data_structure::MotifGraphOP mg) {
    auto motif_names = String();
    for(auto const & n : *mg) {
        motif_names += n->data()->name() + ";";
    }
    return motif_names;
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

    //turn on logging
    base::init_logging();

    //load extra motifs being used
    String base_path = base::base_dir() + "/rnamake/lib/RNAMake/apps/apt_new_interface/resources/";
    resources::Manager::instance().add_motif(base_path+"pRNA_3WJ.pdb", "prna");

    auto app = AptNewInterface();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;

}
