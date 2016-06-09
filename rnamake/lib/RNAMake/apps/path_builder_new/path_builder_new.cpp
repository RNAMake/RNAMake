//
//  path_builder_new.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/29/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "path_builder_new.hpp"


#include "base/file_io.h"
#include "base/settings.h"
#include "util/basic_io.hpp"
#include "util/steric_lookup.hpp"
#include "resources/resource_manager.h"
#include "motif/motif_factory.h"

void
PathBuilderNewApp::setup_options() {
    add_option("pdb", String(""), OptionType::STRING, false);
    add_option("start_bp", String(""), OptionType::STRING, false);
    add_option("end_bp", String(""), OptionType::STRING, false);
    add_option("mg", String(""), OptionType::STRING, false);
    
    add_option("out", String("default.out"), OptionType::STRING, false);
    add_option("no_sterics", 0, OptionType::BOOL, false);
    add_option("write_pdbs", 0, OptionType::BOOL, false);
    add_option("iterate_sterics", 0, OptionType::BOOL, false);
    
    add_cl_options(search_.options(), "search");
    
}

void
PathBuilderNewApp::parse_command_line(
    int argc,
    const char ** argv) {
    
    Application::parse_command_line(argc, argv);
    cl_parser_.assign_options(cl_options_, search_.options(), "search");
    search_.update_var_options();
}

void
PathBuilderNewApp::run() {
    
    if(get_string_option("pdb") != "")     { _setup_from_motif(); }
    else if(get_string_option("mg") != "") { _setup_from_mg(); }
    
    else {
        throw std::runtime_error("not implemented please specify with -pdb");
    }
    
    auto write_pdbs = get_bool_option("write_pdbs");
    if(write_pdbs) {
        //mg_.to_pdb("scaffold.pdb");
        mg_.write_pdbs();
    }
    
    if(get_bool_option("no_sterics")) {
        search_.set_option_value("sterics", false);
    }
    
    
    auto start = mg_.get_node(start_.n_pos)->data()->get_basepair(start_.name)[0]->state();
    auto end = mg_.get_node(end_.n_pos)->data()->get_basepair(end_.name)[0]->state();
    
    
    //auto selector = std::make_shared<MotifStateSelector>();
    //selector->add("unique_twoway");
    //selector->add("ideal_helices_min");
    //selector->connect("unique_twoway", "ideal_helices_min");
    
    //search_.setup(end, start);
    search_.setup(start, end);

    //search_.selector(selector);
    
    int count = 0;
    float dist = 0;
    
    auto beads = Points();
    for(auto & n : mg_) {
        //n->data()->get_beads(n->data()->ends());
        for(auto const & b : n->data()->beads()) {
            if(b.btype() == BeadType::PHOS) { continue; }
            
            dist = b.center().distance(start->d());
            /*if(dist < 40) {
                continue;
            }*/
            beads.push_back(b.center());
            
        }
        
       
    }
    
    auto sl = StericLookup();
    points_to_pdb("beads.pdb", beads);
    mg_.get_node(start_.n_pos)->data()->get_basepair(start_.name)[0]->to_pdb("start.pdb");
    mg_.get_node(end_.n_pos)->data()->get_basepair(end_.name)[0]->to_pdb("end.pdb");
    
    //exit(0);
    sl.add_points(beads);
    
    if(get_bool_option("iterate_sterics")) {
        _iterate_sterics(beads);
        exit(0);
    }
    
    //search_.beads(beads);
    search_.lookup(sl);
    
    std::cout << "num of beads: " << beads.size() << std::endl;
    

    std::ofstream out;
    out.open(get_string_option("out"));
    out << "sol_num\tscore\tsequence\tstructure\ttopology" << std::endl;
    
    int i = -1;
    while(! search_.finished()) {
        i++;
        auto sol = search_.next();
        //std::cout << sol->path().size() << std::endl;
        //exit(0);
        auto mt = sol->to_motif_tree();
        auto ss = mt->secondary_structure();
        
        out << i << "\t" << sol->score() << "\t" << ss->sequence() << "\t" << ss->dot_bracket();
        out << "\t" << mt->topology_to_str() << std::endl;
       
        
        if(write_pdbs) {
            mt->to_pdb("solution."+std::to_string(i) + ".pdb", 1);
        }
        
    }
    
    out.close();

    
}


void
PathBuilderNewApp::_iterate_sterics(Points const & beads) {
    
    bool done = false;
    
    auto start = mg_.get_node(start_.n_pos)->data()->get_basepair(start_.name)[0]->state();
    auto end = mg_.get_node(end_.n_pos)->data()->get_basepair(end_.name)[0]->state();
    
    auto keep_beads = Points();
    auto seen_beads = std::map<int, int>();
    int found = 0;
    
    int count = -1;
    while(! done) {
        count++;
        found = 0;
        auto search = MotifStateSearch();
        search.set_option_value("accept_score", 20);
        search.setup(start, end);
        search.beads(keep_beads);
        
        auto sol = search.next();
        auto mst = sol->to_mst();
        mst->to_motif_tree()->to_pdb("sol."+std::to_string(count) + ".pdb");
        auto sol_beads = Points();
        for(auto const & n : *mst) {
            for(auto const & b : n->data()->cur_state->beads()) {
                sol_beads.push_back(b);
            }
        }
        
        float dist = 0;
        int i = -1;
        for(auto const & b1 : beads) {
            i++;
            if(seen_beads.find(i) != seen_beads.end()) {
                continue;
            }
            for(auto const & b2 : sol_beads) {
                dist = b1.distance(b2);
                if(dist < 2.5) {
                    found = 1;
                    keep_beads.push_back(b1);
                    seen_beads[i] = 1;
                    break;
                }
            }
            
        }
        
        std::cout << keep_beads.size() << std::endl;
        if(found == 0) { break; }
        
        
    }
    
}


void
PathBuilderNewApp::_setup_from_mg() {
    
    auto lines = get_lines_from_file(get_string_option("mg"));
    mg_ =  MotifGraph(lines[0], MotifGraphStringType::MG);
    auto spl = split_str_by_delimiter(lines[1], " ");
    start_ = EndStateInfo{spl[0], std::stoi(spl[1])};
    spl = split_str_by_delimiter(lines[2], " ");
    end_ = EndStateInfo{spl[0], std::stoi(spl[1])};

    
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

    RM::instance().add_motif(get_string_option("pdb"));

    auto m = MotifOP(nullptr);

    try {
        m = get_motif_from_resource_manager(pdb_name);
    }
    
    catch(ResourceManagerException const & e) {
        throw std::runtime_error(
            "cannot load supplied pdb: " + get_string_option("pdb") + " as it is not possible to "
            "build from any basepair ends or it contains no basepair ends. In ability to build "
            "from an end is caused from steric clashes when building a small test helix");
    }
    
    m->to_pdb("motif.pdb");
    
    if(m->ends()[0]->name() == get_string_option("start_bp")) {
        int i = -1;
        int found = 0;
        for(auto const & end : m->ends()) {
            i++;
            if(i == 0) { continue; }
            try {
                m = get_motif_from_resource_manager(pdb_name, "", m->ends()[i]->name());
                found = 1;
                break;
                
            }
            catch(ResourceManagerException const & e) { continue; }
        }
        
        if(!found) {
            throw std::runtime_error(
                "cannot use start_bp: " + start_bp_name + " as its the only end that can be "
                "aligned too. consider swapping start and end basepairs or pick new basepairs " +
                "to build from");
        }
        
    }
    
    auto names = String("");
    for(auto const & end : m->ends()) {
        names += end->name() + " , ";
    }
    
    auto start_bps = m->get_basepair(start_bp_name);
    if(start_bps.size() == 0) {
        throw std::runtime_error(
            start_bp_name + " is not an end of motif given! available ends are: " + names);
    }
    
    auto end_bps = m->get_basepair(end_bp_name);
    if(end_bps.size() == 0) {
        throw std::runtime_error(
            end_bp_name + " is not an end of motif given! available ends are: " + names);
    }
    
    mg_.add_motif(m);
    start_ = EndStateInfo{start_bp_name, 0};
    end_ = EndStateInfo{end_bp_name, 0};
    
   

    
    

    
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