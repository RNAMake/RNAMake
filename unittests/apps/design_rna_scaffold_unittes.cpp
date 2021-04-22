//
// Created by Hassan Abdelsamad on 4/18/21.
//

#include "../command_line_args.hpp"
#include "../common.hpp"

#include "base/cl_option.h"
#include "base/string.h"
#include "base/application.hpp"


class design_rna_scaffold_unittest {

    class DesignRNAScaffold : public base::Application {
    public:

        DesignRNAScaffold();

    public: // application functions

        void
        setup_options() override ;

        void
        parse_command_line(
                int,
                const char **) override ;

        void
        run() override ;
};

TEST_CASE( "Test Application, a class to inherit for excutables") {

    SUBCASE("Test parsing options") {

        auto cla = CommandLineArgs("-test test2 -search.accept_score 5");

        auto app = TestApplication();
        app.parse_command_line(cla.argc, cla.argv());

                CHECK(app.get_string_option("test") == "test2");

                SUBCASE("test parsing options into classes inside application") {
                    CHECK(app.search_.options().get_float("accept_score") == 5);
        }

    }

};