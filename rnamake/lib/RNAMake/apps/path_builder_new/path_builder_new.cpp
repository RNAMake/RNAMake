//
//  path_builder_new.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/29/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "path_builder_new.hpp"


#include "util/file_io.h"
#include "util/settings.h"
#include "util/basic_io.hpp"
#include "util/file_io.h"
#include "resources/resource_manager.h"
#include "motif/motif_factory.h"

void
PathBuilderNewApp::setup_options() {
    add_option("pdb", String(""), OptionType::STRING, false);
    add_option("start_bp", String(""), OptionType::STRING, false);
    add_option("end_bp", String(""), OptionType::STRING, false);
    add_cl_options(optimizer_.options(), "search");

    
    //add_option("f", String(""), OptionType::STRING, true);
    //add_option("break_point", String(""), OptionType::STRING, true);
    //add_option("connections", String(""), OptionType::STRING, false);
    //add_option("o", String("sequences.dat"), OptionType::STRING, false);
}

void
PathBuilderNewApp::parse_command_line(
    int argc,
    const char ** argv) {
    
    Application::parse_command_line(argc, argv);
    cl_parser_.assign_options(cl_options_, optimizer_.options(), "search");
}

void
PathBuilderNewApp::run() {
    
    if(get_string_option("pdb") != "") {
        _setup_from_motif();
    }

}

void
PathBuilderNewApp::_setup_from_motif() {
    
    auto pdb_name = filename(get_string_option("pdb"));
    auto start_bp_name = get_string_option("start_bp");
    auto end_bp_name = get_string_option("end_bp");

    if(start_bp_name == "" || end_bp_name == "") {
        throw std::runtime_error(
            "must supply the name of the start_bp and end_bp using option -start_bp and "
            "-end_bp respectively when using -pdb option");
    }
    
    ResourceManager::getInstance().add_motif(get_string_option("pdb"));
    auto m = ResourceManager::getInstance().get_motif(pdb_name);
    m->to_pdb("motif.pdb");
    
        //if(m->ends()[0]->name() == get_string_option("start_bp")) {
    //    m = ResourceManager::getInstance().get_motif(pdb_name, "", m->ends()[3]->name());
    //}
    
    
    
    /*auto beads = m->get_beads(m->ends());
    auto centers = Points();
    
    for(auto const & b : beads) {
        if(b.btype() == BeadType::PHOS) { continue; }
        centers.push_back(b.center());
    }
    points_to_pdb("points.pdb", centers);
    */
    
    auto start = m->get_basepair(get_string_option("start_bp"))[0]->state();
    auto end   = m->get_basepair(get_string_option("end_bp"))[0]->state();

    auto search = MotifStateSearch();
    search.setup(start, end);
    //search.beads(centers);
    search.set_option_value("verbose", true);
    
    auto sol = search.next();
    sol->to_pdb("sol.pdb", 1);
    
    

    
}

int main(int argc, const char * argv[]) {
    
    auto app = PathBuilderNewApp();
    try {
        app.setup_options();
        app.parse_command_line(argc, argv);
        app.run();
    }
    catch(std::runtime_error const & e) {
        std::cout << e.what() << std::endl;
    }
    return 0;
    
}