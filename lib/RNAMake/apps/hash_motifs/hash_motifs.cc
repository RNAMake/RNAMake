//
// Created by Joseph Yesselman on 3/15/19.
//


#include "base/backtrace.hpp"
#include "base/log.h"
#include "base/settings.h"
#include "math/euler.h"
#include "hash_motifs/hash_motifs.h"
#include <util/cartesian_product.h>
#include <motif/motif_state_aligner.h>
#include <motif_search/monte_carlo/scorer.h>
#include <motif_search/exhaustive/search.h>

HashMotifs::HashMotifs() { }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// app functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
HashMotifs::setup_options() {

}

void
HashMotifs::parse_command_line(
        int argc,
        const char ** argv) {

}

void
HashMotifs::run() {
    auto & rm = resources::Manager::instance();

    String base_path = base::base_dir() + "/rnamake/lib/RNAMake/apps/simulate_tectos/resources/";
    rm.add_motif(base_path+"GAAA_tetraloop");

    //_old_buildup();
    //exit(0);

    auto ttr = rm.motif("GAAA_tetraloop", "", "A229-A245");
    auto start_pos = ttr->get_end_index("A222-A251");
    auto end_pos = ttr->get_end_index("A149-A154");

    auto start_bp = ttr->ends()[end_pos]->state();
    auto end_bp = ttr->ends()[start_pos]->state();

    auto aligner = motif::MotifStateAligner();
    auto h = rm.motif_state("HELIX.IDEAL.2");
    aligner.get_aligned_motif_state(start_bp, h);

    auto sol_template = motif_search::SolutionTopologyTemplate();
    sol_template.add_library("twoway");
    sol_template.add_library("flex_helices", data_structure::NodeIndexandEdge{0, 1});
    sol_template.add_library("twoway", data_structure::NodeIndexandEdge{1, 1});
    sol_template.add_library("flex_helices", data_structure::NodeIndexandEdge{2, 1});

    auto factory = motif_search::SolutionToplogyFactory();
    auto sol_toplogy = factory.generate_toplogy(sol_template);
    auto filter = std::make_shared<motif_search::NoExclusionFilter>();

    auto scorer = std::make_shared<motif_search::exhaustive::DefaultScorer>();
    auto search = motif_search::exhaustive::Search(scorer, *sol_toplogy, filter);

    auto lookup = std::make_shared<util::StericLookupNew>();
    auto p = std::make_shared<motif_search::Problem>(h->end_states()[1], end_bp, lookup, false);

    search.setup(p);
    auto sol = search.next();


    //generate_pairwise_map();
    //exit(0);

    /*auto bb = math::BoundingBox(math::Point(-100, -100, -100), math::Point(100, 100, 100));
    auto bin_widths = math::Real6{0.5, 0.5, 0.5, 10, 10, 10};
    auto hash = MotifPathHash(bb, bin_widths);
    */

    /*exit(0);

    auto values = math::Real6();
    auto euler = math::Vector();
    int i = -1, j = -1;
    for(auto const & ms1 : flex_helices_lib) {
        i++;
        j = -1;
        for(auto const & ms2 : twoway_lib ) {
            j++;
            auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
            msg->set_option_value("sterics", false);
            msg->add_state(ms1);
            msg->add_state(ms2);

            auto & t = msg->get_node(1)->data()->get_end_state(1)->d();
            auto r = msg->get_node(1)->data()->get_end_state(1)->r();
            math::calc_euler(r, euler);

            for(int i = 0; i < 3; i++) {
                euler[i] = euler[i]*180/M_PI;
                if(euler[i] > 180) {
                    euler[i] -= 360;
                }
            }

            values[0] = t[0];
            values[1] = t[1];
            values[2] = t[2];
            values[3] = euler[0];
            values[4] = euler[1];
            values[5] = euler[2];

            hash.add(values, Ints{i, j});
        }
    }

    i = -1, j = -1;
    for(auto const & ms1 : flex_helices_lib) {
        i++;
        j = -1;
        for(auto const & ms2 : twoway_lib ) {
            j++;
            auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
            msg->set_option_value("sterics", false);
            msg->add_state(ms1);
            msg->add_state(ms2);

            auto & t = msg->get_node(1)->data()->get_end_state(1)->d();
            auto r = msg->get_node(1)->data()->get_end_state(1)->r();
            math::calc_euler(r, euler);

            for(int i = 0; i < 3; i++) {
                euler[i] = euler[i]*180/M_PI;
                if(euler[i] > 180) {
                    euler[i] -= 360;
                }
            }

            values[0] = t[0];
            values[1] = t[1];
            values[2] = t[2];
            values[3] = euler[0];
            values[4] = euler[1];
            values[5] = euler[2];

            auto results = hash.search(values);
            for(auto const & result : results) {
                for(auto const & e : result) {
                    std::cout << e << " ";
                }
                std::cout << std::endl;
            }
            exit(0);
        }
    }

    exit(0);

    std::ofstream out;
    out.open("test.out", std::ios::binary);
    hash.output_binary(out);
    out.close();

    std::ifstream in;
    in.open("test.out", std::ios::binary);
    auto hash_2 = MotifPathHash(in);*/

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
HashMotifs::_generate_pairwise_map() {

    auto flex_helices_lib =  resources::MotifStateSqliteLibrary("flex_helices");
    auto twoway_lib =  resources::MotifStateSqliteLibrary("twoway");

    flex_helices_lib.load_all();
    twoway_lib.load_all();


    auto values = math::Real6();
    auto euler = math::Vector();
    auto final_bps = structure::BasepairStateOPs();
    int i = -1, j = -1;
    for(auto const & ms2 : twoway_lib ) {
        j++;
        auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
        msg->set_option_value("sterics", false);
        msg->add_state(ms2);

        auto & t = msg->get_node(0)->data()->get_end_state(1)->d();
        auto r = msg->get_node(0)->data()->get_end_state(1)->r();

        final_bps.push_back(std::make_shared<structure::BasepairState>(*msg->get_node(0)->data()->get_end_state(1)));
    }

    i = -1; j = -1;
    auto map = std::vector<Floats>(final_bps.size());
    for(auto const & bp1 : final_bps) {
        i++;
        map[i] = Floats(final_bps.size());
        for(auto const & bp2 : final_bps) {
            j++;
            map[i][i] = bp1->diff(bp2);
        }
    }
}

void
HashMotifs::_old_buildup() {
    auto & rm = resources::Manager::instance();
    auto flex_helices_lib =  resources::MotifStateSqliteLibrary("flex_helices");
    auto twoway_lib =  resources::MotifStateSqliteLibrary("twoway");

    flex_helices_lib.load_all();
    twoway_lib.load_all();

    auto flex_helices_ms_1 = motif::MotifStateOPs();
    auto flex_helices_ms_2 = motif::MotifStateOPs();
    auto flex_helices_ms_3 = motif::MotifStateOPs();
    for(auto const & ms : flex_helices_lib) {
        if(ms->size() > 14) { continue; }
        flex_helices_ms_1.push_back(ms);
        flex_helices_ms_2.push_back(std::make_shared<motif::MotifState>(*ms));
        flex_helices_ms_3.push_back(std::make_shared<motif::MotifState>(*ms));
    }

    auto twoway_ms_1 = motif::MotifStateOPs();
    auto twoway_ms_2 = motif::MotifStateOPs();
    for(auto const & ms : twoway_lib) {
        twoway_ms_1.push_back(ms);
        twoway_ms_2.push_back(std::make_shared<motif::MotifState>(*ms));
    }

    auto scorer = motif_search::monte_carlo::GreedyScorer();
    auto ttr = resources::Manager::instance().motif("GAAA_tetraloop", "", "A229-A245");
    //ttr->to_pdb("test.0.pdb", 1, 1);

    auto start_pos = ttr->get_end_index("A222-A251");
    auto end_pos = ttr->get_end_index("A149-A154");

    auto start_bp = ttr->ends()[end_pos]->state();
    auto end_bp = ttr->ends()[start_pos]->state();

    auto aligner = motif::MotifStateAligner();
    auto h = rm.motif_state("HELIX.IDEAL.2");
    aligner.get_aligned_motif_state(start_bp, h);
    //rm.get_motif_from_state(h)->to_pdb("test.1.pdb", 1, 1);


    scorer.set_target(end_bp, false);

    auto best = 10000.0f;
    auto dist = 0;
    auto score = 0.0f;
    auto count = 0;
    for(auto const & ms1 : twoway_ms_1) {
        std::cout << ms1->name() << std::endl;
        aligner.get_aligned_motif_state(h->end_states()[1], ms1);
        for(auto const & ms2 : flex_helices_ms_1) {
            aligner.get_aligned_motif_state(ms1->end_states()[1], ms2);
            for(auto const & ms3 : twoway_ms_2) {
                aligner.get_aligned_motif_state(ms2->end_states()[1], ms3);
                for(auto const & ms4 : flex_helices_ms_2) {
                    aligner.get_aligned_motif_state(ms3->end_states()[1], ms4);
                    /*for(auto const & ms5 : flex_helices_ms_2) {
                        aligner.get_aligned_motif_state(ms4->end_states()[1], ms5);
                        score = scorer.score(*ms5->end_states()[1]);
                    }*/
                    score = scorer.score(*ms4->end_states()[1]);
                    count += 1;
                    //if(count % 1000000 == 0) { std::cout << count << std::endl; }

                    //std::cout << ms4->end_states()[1]->d() << std::endl;
                    //std::cout << score << std::endl;
                    /*for(auto const & ms : motif::MotifStateOPs{ms1, ms2, ms3, ms4}) {
                        std::cout << ms->name() << " ";
                    }
                    std::cout << std::endl;
                    */

                    if(score < best) {
                        best = score;
                        //std::cout << best << std::endl;
                        /*rm.get_motif_from_state(ms1)->to_pdb("test.2.pdb", 1, 1);
                        rm.get_motif_from_state(ms2)->to_pdb("test.3.pdb", 1, 1);
                        rm.get_motif_from_state(ms3)->to_pdb("test.4.pdb", 1, 1);
                        rm.get_motif_from_state(ms4)->to_pdb("test.5.pdb", 1, 1);*/

                    }

                    if(score < 10) {
                        exit(0);
                        /*auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
                        msg->set_option_value("sterics", false);
                        msg->add_state(ttr->get_state());
                        msg->add_state(h);
                        msg->add_state(ms1);
                        msg->add_state(ms2);
                        msg->add_state(ms3);
                        msg->add_state(ms4);
                        msg->add_connection(0, 5, "A222-A251", ms4->end_names()[1]);
                        msg->to_motif_graph()->to_pdb("design."+std::to_string(count)+".pdb", 1, 1);
                        count += 1;
                        std::cout << score << std::endl;*/
                    }

                }
            }
        }
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

    auto app = HashMotifs();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;

}