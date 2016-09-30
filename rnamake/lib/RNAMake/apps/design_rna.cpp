//
//  design_rna.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/26/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "design_rna.hpp"

#include "base/backtrace.hpp"
#include "resources/resource_manager.h"
#include "motif_tools/segmenter.h"



void
DesignRNAApp::setup_options() {
    // start from pdb
    add_option("pdb", String(""), OptionType::STRING, false);
    add_option("start_bp", String(""), OptionType::STRING, false);
    add_option("end_bp", String(""), OptionType::STRING, false);

    add_cl_options(search_.options(), "search");
}

void
DesignRNAApp::parse_command_line(
    int argc,
    const char ** argv) {
    
    Application::parse_command_line(argc, argv);
    cl_parser_.assign_options(cl_options_, search_.options(), "search");
    search_.update_var_options();

}


void
DesignRNAApp::_setup_sterics() {
    auto beads = Points();
    for(auto & n : mg_) {
        for(auto const & b : n->data()->beads()) {
            if(b.btype() == BeadType::PHOS) { continue; }
            beads.push_back(b.center());
            
        }
        
        for(auto const & b : n->data()->protein_beads()) { beads.push_back(b.center()); }
    }
    lookup_.add_points(beads);
}


void
DesignRNAApp::_setup_from_pdb() {
    auto struc = RM::instance().get_structure(get_string_option("pdb"), "scaffold");
    
    std::cout << "DESIGN RNA: loaded pdb from file: " << get_string_option("pdb") << std::endl;

    auto start_bp_name = get_string_option("start_bp");
    auto end_bp_name = get_string_option("end_bp");
    
    if(start_bp_name == "" || end_bp_name == "") {
        throw DesignRNAAppException(
            "must supply the name of the start_bp and end_bp using option -start_bp and "
            "-end_bp respectively when using -pdb option");
    }
    
    auto start_bps = struc->get_basepair(start_bp_name);
    auto end_bps = struc->get_basepair(end_bp_name);
    
    if(start_bps.size() == 0) {
        throw DesignRNAAppException("cannot find start basepair: " + start_bp_name);
    }
    
    if(start_bps.size() > 1) {
        throw DesignRNAAppException(
            "start basepair name: " + start_bp_name + " is not unique name, multiple basepairs "
            "have this name please pick another basepair or reformat your pdb");
    }
    
    if(end_bps.size() == 0) {
        throw DesignRNAAppException("cannot find end basepair: " + end_bp_name);
    }
    
    if(end_bps.size() > 1) {
        throw DesignRNAAppException(
            "end basepair name: " + end_bp_name + " is not unique name, multiple basepairs "
            "have this name please pick another basepair or reformat your pdb");
    }

    if(start_bps[0]->bp_type() != "cW-W") {
        throw DesignRNAAppException(
            "start basepair is not a watson and crick basepair building from non-canonical "
            "leads to bizzare effects");
    }
    
    if(end_bps[0]->bp_type() != "cW-W") {
        throw DesignRNAAppException(
            "ends basepair is not a watson and crick basepair building from non-canonical "
            "leads to bizzare effects");
    }
    
    
    auto bps = BasepairOPs{start_bps[0], end_bps[0]};
    
    auto segmenter = Segmenter();
    auto segments = segmenter.apply(struc, bps);
    
    std::cout << "DESIGN RNA: segmentation success" << std::endl;
    std::cout << "DESIGN RNA: removed rna segment=removed.pdb" << std::endl;
    std::cout << "DESIGN RNA: remaining rna segment=remaining.pdb" << std::endl;

    segments->remaining->to_pdb("remaining.pdb");
    segments->removed->to_pdb("removed.pdb");
    
    //auto new_struc = std::make_shared<RNAStructure>(*struc);
    segments->remaining->block_end_add(-1);
    mg_.add_motif(segments->remaining);
    start_ = EndStateInfo{start_bp_name, 0};
    end_ = EndStateInfo{end_bp_name, 0};
}

void
DesignRNAApp::run() {
    
    if(get_string_option("pdb") != "") { _setup_from_pdb(); }
    
    auto start = mg_.get_node(start_.n_pos)->data()->get_basepair(start_.name)[0]->state();
    auto end = mg_.get_node(end_.n_pos)->data()->get_basepair(end_.name)[0]->state();
    _setup_sterics();
    
    search_.setup(start, end);
    search_.lookup(lookup_);
    
    auto sols = std::vector<MotifTreeOP>();
    
    int i = -1;
    while(! search_.finished()) {
        i++;
        auto sol = search_.next();
        
        auto mt = sol->to_motif_tree();
        sols.push_back(mt);
        
    }
    
    
    
}


int main(int argc, const char * argv[]) {
    std::set_terminate(print_backtrace);
    
    auto app = DesignRNAApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();
    return 0;
    
}