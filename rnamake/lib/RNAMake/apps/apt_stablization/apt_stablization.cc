//
// Created by Joseph Yesselman on 4/12/18.
//

#include "base/backtrace.hpp"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_state_graph.hpp"
#include "motif_state_search/motif_state_monte_carlo.h"
#include "sequence_optimizer/sequence_optimizer_3d.hpp"
#include "apt_stablization/apt_stablization.h"



APTStablization::APTStablization() : Application() {}

void
APTStablization::setup_options() {
    add_option("pdb", String(""), OptionType::STRING, true);

    // general options
    add_option("out_file", "default.out", OptionType::STRING, false);
    add_option("score_file", "default.scores", OptionType::STRING, false);
    add_option("designs", 1, OptionType::INT, false);

}

void
APTStablization::parse_command_line(
        int argc,
        const char **argv) {

    Application::parse_command_line(argc, argv);
}


void
APTStablization::run() {

    // add motif to resource manager
    RM::instance().add_motif(get_string_option("pdb"), "aptamer", MotifType::TWOWAY);
    std::cout << "APT STABLIZATION: loaded pdb from file: " << get_string_option("pdb") << std::endl;

    auto m = RM::instance().motif("aptamer");
    if(m->ends().size() != 2) {
        throw APTStablizationException("aptamer does not have not have two basepair ends");
    }

    //auto start_path = "flex_helices,twoway,flex_helices,twoway,flex_helices,aptamer,flex_helices";
    auto start_path = "flex_helices,aptamer,flex_helices,twoway,flex_helices,twoway,flex_helices";
    auto ms_libraries = _get_libraries(start_path);

    // setup initial graph
    auto ttr = RM::instance().motif("GAAA_tetraloop", "", "A229-A245");
    auto bp_step_1 = RM::instance().bp_step("GC_LL_GC_RR");
    auto bp_step_2 = RM::instance().bp_step("CG_LL_CG_RR");
    auto msg = std::make_shared<MotifStateGraph>();
    msg->add_state(bp_step_1->get_state());
    msg->add_state(bp_step_2->get_state());
    msg->add_state(ttr->get_state());

    _setup_sterics(msg);

    auto start_end_pos = msg->get_node(2)->data()->get_end_index("A222-A251");
    auto end_end_pos = msg->get_node(2)->data()->get_end_index("A149-A154");

    bool target_an_aligned_end = false;
    if(end_end_pos == msg->get_node(2)->data()->block_end_add()) {
        target_an_aligned_end = true;
    }

    auto mc = MotifStateMonteCarlo(ms_libraries);
    mc.setup(msg, 2, 2, start_end_pos, end_end_pos, target_an_aligned_end);
    mc.lookup(lookup_);
    mc.start();
    mc.set_option_value("max_solutions", 10000);

    auto solutions = 0;

    while(! mc.finished()) {
        auto mg = mc.next();
        if(mg == nullptr) {
            std::cout << "no solutions" << std::endl;
            exit(0);
        }
        //mg->to_pdb("design." + std::to_string(0) + ".pdb", 1, 1, 1);
        for(auto & n : *mg) {
            if(n->data()->name().substr(0,5) == "HELIX") {
                n->data()->mtype(MotifType::HELIX);
            }
            std::cout << n->data()->name() << " " << n->data()->mtype() << std::endl;
        }
        mg->to_pdb("test.pdb", 1);
        mg->replace_ideal_helices();

        auto end_node = mg->get_node(2);
        auto target_state = end_node->data()->ends()[end_end_pos]->state();
        auto end_i = end_end_pos;
        auto partner = end_node->connections()[end_i]->partner(end_node->index());
        std::cout << partner->index() << std::endl;
        //auto scorer = std::make_shared<ExternalTargetScorer>(target_state, partner->index(), 1, target_an_aligned_end);
        auto scorer = std::make_shared<InternalTargetScorer>(2, end_end_pos, partner->index(), 1, target_an_aligned_end);

        auto optimizer = SequenceOptimizer3D();
        optimizer.set_option_value("verbose", true);
        optimizer.set_option_value("cutoff", 5.0f);

        auto sols = optimizer.get_optimized_sequences(mg, scorer);
        if(sols[0]->dist_score > 5) {
            continue;
        }

        auto copy_mg = std::make_shared<MotifGraph>(*mg);
        auto dss = copy_mg->designable_secondary_structure();
        dss->replace_sequence(sols[0]->sequence);
        copy_mg->replace_helical_sequence(dss);
        copy_mg->to_pdb("design." + std::to_string(solutions) + ".pdb", 1, 1, 1);
        exit(0);



    }

    //mc.run();



}


std::vector<MotifStateOPs>
APTStablization::_get_libraries(
        String const & motif_path) {
    auto spl = split_str_by_delimiter(motif_path, ",");
    auto i = 0;
    auto libraries = std::vector<MotifStateOPs>();
    auto motif_states = MotifStateOPs();
    for(auto const & name : spl) {
        if(name.length() < 2) { continue; }
        if(name == "ideal_helices_min" || name == "unique_twoway" || name == "tcontact" ||
           name == "twoway" || name == "flex_helices") {
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
                auto m_new = RM::instance().motif(name, "", end->name());
                m_new->new_res_uuids();
                //m_new->ends()[1]->flip();
                //std::cout << m_new->name() << std::endl;
                motif_states.push_back(m_new->get_state());
            }
            libraries.push_back(motif_states);
        }
    }
    return libraries;
}

void
APTStablization::_setup_sterics(
        MotifStateGraphOP msg) {
    auto beads = Points();
    for (auto & n : *msg) {
        for(auto const & b : n->data()->cur_state->beads()) {
            beads.push_back(b);
        }
    }
    lookup_.add_points(beads);
}

int main(int argc, const char *argv[]) {
    //must add this for all apps!
    std::set_terminate(print_backtrace);

    String base_path = base_dir() + "/rnamake/lib/RNAMake/apps/simulate_tectos/resources/";
    RM::instance().add_motif(base_path+"GAAA_tetraloop");

    auto app = APTStablization();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;

}
