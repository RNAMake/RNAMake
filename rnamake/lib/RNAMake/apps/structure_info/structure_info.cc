//
// Created by Joseph Yesselman on 11/21/18.
//

#include "base/backtrace.hpp"
#include "resources/resource_manager.h"
#include "structure_info/structure_info.h"

StructureInfoApp::StructureInfoApp() : Application() {

}

// application setups functions ////////////////////////////////////////////////////////////////////

void
StructureInfoApp::setup_options() {
    add_option("pdb", "", OptionType::STRING, true);
    add_option("basepairs", false, OptionType::BOOL, false);
}

void
StructureInfoApp::parse_command_line(
        int argc,
        const char ** argv) {
    Application::parse_command_line(argc, argv);
}

void
StructureInfoApp::run() {
    auto struc = RM::instance().get_structure(get_string_option("pdb"), "test");

    std::cout << "loaded pdb: " << get_string_option("pdb") << std::endl;
    std::cout << "sequence: ";
    size_t i = 0;
    for(auto const & c : struc->chains()) {
        i++;
        for(auto const & r : c->residues()) {
            std::cout << r->name();
        }
        if(i != struc->chains().size()) {  std::cout << "&"; }
    }
    std::cout << std::endl;

    std::cout << "basepairs: " << std::endl;
    for(auto const & bp : struc->basepairs()) {
        std::cout << bp->name() << " " << bp->bp_type() << std::endl;
    }


}

int main(int argc, const char * argv[]) {
    //must add this for all apps!
    std::set_terminate(print_backtrace);

    auto app = StructureInfoApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;
}