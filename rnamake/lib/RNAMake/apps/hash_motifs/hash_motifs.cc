//
// Created by Joseph Yesselman on 3/15/19.
//


#include "base/backtrace.hpp"
#include "base/log.h"
#include "base/settings.h"
#include "math/euler.h"
#include "hash_motifs/hash_motifs.h"

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
    auto bb = math::BoundingBox(math::Point(-100, -100, -100), math::Point(100, 100, 100));
    auto bin_widths = math::Real6{0.5, 0.5, 0.5, 10, 10, 10};
    auto hash = MotifPathHash(bb, bin_widths);

    auto flex_helices_lib =  resources::MotifStateSqliteLibrary("flex_helices");
    auto twoway_lib =  resources::MotifStateSqliteLibrary("twoway");

    flex_helices_lib.load_all(10);
    twoway_lib.load_all(10);


    auto values = math::Real6();
    auto euler = math::Vector();

    for(auto const & ms1 : flex_helices_lib) {
        for(auto const & ms2 : twoway_lib ) {
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

            hash.add_value(values);
        }
    }

    std::ofstream out;
    out.open("test.out", std::ios::binary);
    hash.output_binary(out);
    out.close();

    std::ifstream in;
    in.open("test.out", std::ios::binary);
    auto hash_2 = MotifPathHash(in);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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