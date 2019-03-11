//
// Created by Joseph Yesselman on 2/23/19.
//

#include "base/backtrace.hpp"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_state_graph.hpp"
#include "motif_state_search/motif_state_monte_carlo.h"
#include "sequence_optimizer/sequence_optimizer_3d.hpp"
#include "apt_new_interface/apt_new_interface.h"


AptNewInterface::AptNewInterface() {}

void
AptNewInterface::setup_options() {
    add_option("scaffold", String(""), base::OptionType::STRING, true);
    add_option("docked_motif", String(""), base::OptionType::STRING, true);
    add_option("scaffold_end", String(""), base::OptionType::STRING, true);

    // general options
    add_option("out_file", "default.out", base::OptionType::STRING, false);
    add_option("score_file", "default.scores", base::OptionType::STRING, false);
    add_option("designs", 1, base::OptionType::INT, false);

}

void
AptNewInterface::parse_command_line(
        int argc,
        const char **argv) {

    base::Application::parse_command_line(argc, argv);
}

void
AptNewInterface::run() {
    auto mf = MotifFactory();

    // add new motifs to resource manager
    auto scaffold_rm = RM::instance().get_structure(get_string_option("scaffold"), "scaffold", 3);
    auto scaffold_m = std::make_shared<Motif>(*scaffold_rm);
    scaffold_m->name("scaffold");
    scaffold_m->mtype(MotifType::TCONTACT);
    mf._setup_secondary_structure(scaffold_m);
    RM::instance().register_motif(scaffold_m);

    //RM::instance().add_motif(get_string_option("scaffold"), "scaffold", MotifType::TCONTACT);
    auto rs = RM::instance().get_structure(get_string_option("docked_motif"), "docked_motif");

    auto prna = RM::instance().motif("prna", "", "A7-C10");
    auto scaffold = RM::instance().motif("scaffold", "", get_string_option("scaffold_end"));
    auto docked_motif = std::make_shared<Motif>(*rs);
    docked_motif->mtype(MotifType::HAIRPIN);
    mf._setup_secondary_structure(docked_motif);

    RM::instance().register_motif(docked_motif);

    auto msg = std::make_shared<MotifStateGraph>();
    msg->set_option_value("sterics", false);
    msg->add_state(scaffold->get_state());
    msg->add_state(RM::instance().motif_state("HELIX.IDEAL.1"), 0, 2);
    msg->add_state(prna->get_state());
    msg->add_state(docked_motif->get_state(), -1, -1, 1);
    msg->increase_level();

    auto msg_copy = std::make_shared<MotifStateGraph>(*msg);

    lookup_ = StericLookup(1.0, 5.0, 7);
    _setup_sterics(msg_copy);

    auto start_path_1 = String("flex_helices,twoway,flex_helices");
    auto start_path_2 = String("flex_helices,twoway,flex_helices,twoway,flex_helices");

    //auto start_path_1 = String("ideal_helices_min,twoway,ideal_helices_min");
    //auto start_path_2 = String("ideal_helices_min,twoway,ideal_helices_min,twoway,ideal_helices_min");

    auto ms_libraries_1 = _get_libraries(start_path_1);
    auto ms_libraries_2 = _get_libraries(start_path_2);

    auto start_end_pos = msg->get_node(2)->data()->get_end_index("B14-C7");
    auto end_end_pos = msg->get_node(0)->data()->get_end_index("A41-A87");

    auto mc_1 = MotifStateMonteCarlo(ms_libraries_1);
    mc_1.setup(msg_copy, 2, 0, start_end_pos, end_end_pos, false);
    mc_1.set_option_value("accept_score", 10);
    mc_1.set_option_value("max_solutions", 1);
    mc_1.lookup(lookup_);
    mc_1.start();

    auto sol_1 = mc_1.next_state();
    if(sol_1 == nullptr) {
        std::cout << " no solution 1" << std::endl;
        exit(0);
    }

    else {
        std::cout << " first solution found! " << std::endl;
    }

    auto sol_1_msg = sol_1->msg;
    /*for(auto it = sol_1_msg->node_begin();
             it != sol_1_msg->node_end();
             it++){

        auto n = *(it);
        std::cout << n->index() << " : ";
        for(auto const & c : n->connections()) {
            if(c == nullptr) { continue; }
            std::cout << c->partner(n->index())->index() << " ";
        }
        std::cout << std::endl;
    }*/

    auto next_msg = 0;

    auto beads = math::Points();
    for (auto & n : *msg) {
        for(auto const & b : n->data()->cur_state->beads()) {
            beads.push_back(b);
        }
    }
    lookup_.add_points(beads);

    start_end_pos = msg->get_node(2)->data()->get_end_index("A10-B9");
    end_end_pos = msg->get_node(3)->data()->get_end_index("A60-A65");

    auto mc_2 = MotifStateMonteCarlo(ms_libraries_2);
    mc_2.setup(sol_1_msg, 2, 3, start_end_pos, end_end_pos, true);
    mc_2.set_option_value("accept_score", 15);
    mc_2.set_option_value("max_solutions", 1);
    mc_2.lookup(lookup_);
    mc_2.start();

    auto sol_2 = mc_2.next();
    if(sol_2 == nullptr) {
        std::cout << " no solution 2" << std::endl;
        exit(0);
    }

    else {
        std::cout << " second solution found! " << std::endl;
    }


    sol_2->mg->replace_ideal_helices();
    sol_2->mg->write_pdbs();
    sol_2->mg->to_pdb("test.pdb", 1, 1);
    auto p = sol_2->mg->secondary_structure();
    std::cout << p->sequence() << std::endl;
    std::cout << p->dot_bracket() << std::endl;

    for(auto const & n : *sol_2->mg) {
        std::cout << n->index() << " : ";
        for(auto const & c : n->connections()) {
            if(c == nullptr) { continue; }
            std::cout << c->partner(n->index())->index() << " ";
        }
        std::cout << std::endl;
    }

    std::ofstream out;
    out.open("multi_path_solution.out");
    out << sol_2->mg->to_str();
    out.close();




}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<MotifStateOPs>
AptNewInterface::_get_libraries(
        String const & motif_path) {
    auto spl = base::split_str_by_delimiter(motif_path, ",");
    auto i = 0;
    auto libraries = std::vector<MotifStateOPs>();
    auto motif_states = MotifStateOPs();
    for(auto const & name : spl) {
        if(name.length() < 2) { continue; }
        if(name == "ideal_helices_min" || name == "unique_twoway" || name == "tcontact" ||
           name == "twoway" || name == "flex_helices" || name == "existing") {
            auto ms_lib =  MotifStateSqliteLibrary(name);
            ms_lib.load_all();
            motif_states = MotifStateOPs();
            if(name == "flex_helices") {
                for (auto const & ms : ms_lib) {
                    if(ms->size() < 20) {
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
                throw AptNewInterfaceException("no viable aptamer conformations");
            }

            libraries.push_back(motif_states);
        }
    }
    return libraries;
}

void
AptNewInterface::_setup_sterics(
        MotifStateGraphOP msg) {
    auto beads = math::Points();
    for (auto & n : *msg) {
        for(auto const & b : n->data()->cur_state->beads()) {
            beads.push_back(b);
        }
    }
    lookup_.add_points(beads);
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
    String base_path = base::base_dir() + "/rnamake/lib/RNAMake/apps/apt_new_interface/resources/";
    RM::instance().add_motif(base_path+"pRNA_3WJ.pdb", "prna");

    auto app = AptNewInterface();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;

}