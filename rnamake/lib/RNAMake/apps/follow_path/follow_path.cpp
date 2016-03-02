//
//  follow_path.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/21/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "follow_path.h"
#include "util/file_io.h"
#include "util/basic_io.hpp"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_topology.h"
#include "motif_data_structures/motif_graph.h"
#include "motif_data_structures/motif_state_tree.h"
#include "motif_state_search/motif_state_search.h"

void
PathBuilder::setup_options() {
    options_.add_option("path", String(""), OptionType::STRING);
    options_.add_option("mg", String(""), OptionType::STRING);
    options_.lock_option_adding();
    update_var_options();
}

void
PathBuilder::setup(
    CommandLineOptions const & cmd_opts) {
    
    for(auto const & opt: cmd_opts) {
        if(! has_option(opt->name())) { continue; }
        
        if     (opt->type() == OptionType::INT) {
            set_option_value(opt->name(), opt->get_int());
        }
        else if(opt->type() == OptionType::FLOAT) {
            set_option_value(opt->name(), opt->get_float());
        }
        else if(opt->type() == OptionType::BOOL) {
            set_option_value(opt->name(), opt->get_bool());
        }
        else if(opt->type() == OptionType::STRING) {
            set_option_value(opt->name(), opt->get_string());
        }
    }
    
}

void
PathBuilder::build() {
    auto lines = get_lines_from_file(get_string_option("path"));
    auto path_points = vectors_from_str(lines[0]);
    
    auto mg = MotifGraphOP(nullptr);
    int n1;
    String end_name = "";
    if(get_string_option("mg").size() > 0) {
        auto lines2 = get_lines_from_file(get_string_option("mg"));
        mg = std::make_shared<MotifGraph>(lines2[0]);
        
        auto spl = split_str_by_delimiter(lines2[1], " ");
        n1 = std::stoi(spl[0]);
        end_name = spl[1];
    }
    
    auto segments = _get_segments(path_points);
    auto pathes = _get_sub_pathes(segments);
    
    auto mt = graph_to_tree(mg, mg->get_node(0));
    auto mst = std::make_shared<MotifStateTree>(mt);
    
    auto pf = PathFollower();
    pf.setup(pathes[0], mg, n1, end_name);
    
    auto sol = pf.solutions();
    mst->add_mst(sol[0]->to_mst(), 0, -1, end_name);

    pf = PathFollower();
    
    
}


std::vector<Points>
PathBuilder::_get_segments(
    Points const & path) {
    
    auto direction = Vector();
    auto new_direction = Vector();
    auto segment = Points();
    auto segments = std::vector<Points>();
    
    int i = -1;
    for(auto const & p : path) {
        i++;
        if(i < 2) {
            segment.push_back(p);
            if(i == 1) {
                direction = segment[0] - segment[1];
                direction = direction.normalize();
            }
            continue;
        }
        
        new_direction = segment.back() - p;
        new_direction = new_direction.normalize();
        
        if(direction.distance(new_direction) > 0.1) {
            segments.push_back(segment);
            segment = Points();
            segment.push_back(p);
            direction = new_direction;
        }
        else {
            segment.push_back(p);
        }
    }
    
    if(segment.size() > 0) {
        segments.push_back(segment);
    }
    
    return segments;
    
}


std::vector<Points>
PathBuilder::_get_sub_pathes(
    std::vector<Points> const & segments) {
    
    auto pathes = std::vector<Points>();
    for(int i = 1; i < segments.size(); i++) {
        auto p = Points();
        if(i == 1) {
            for(auto const & e : segments[0]) { p.push_back(e); }
        }
        else {
            auto half = (int)segments[i-1].size() / 2;
            for(int j = half-1; j < segments[i-1].size(); j++) { p.push_back(segments[i-1][j]); }
        }
        
        if(i == segments.size()-1) {
            for(auto const & e : segments[i]) { p.push_back(e); }
        }
        else {
            auto half = (int)segments[i].size() / 2;
            for(int j = 0; j < half; j++) { p.push_back(segments[i][j]); }
        }
        
        pathes.push_back(p);
    }
    
    return pathes;
}

//MAIN
/**************************************************************************************************/


CommandLineOptions
parse_command_line(
    int argc,
    const char ** argv) {
    
    CommandLineOptions cl_opts;
    auto search = MotifStateSearch();
    cl_opts.add_options(search.options());
    cl_opts.add_option("path", String(""), OptionType::STRING, false);
    cl_opts.add_option("mg", String(""), OptionType::STRING, false);
    cl_opts.add_option("only_one", 1, OptionType::INT, false);
    cl_opts.add_option("full_path", false, OptionType::BOOL, false);
    cl_opts.parse_command_line(argc, argv);
    
    return cl_opts;
    
}


int main(int argc, const char * argv[]) {
    auto cmd_opts = parse_command_line(argc, argv);

    //load TTR
    auto base_dir = String("/Users/josephyesselman/projects/RNAMake/rnamake/lib/RNAMake/apps/mini_ttr");
    ResourceManager::getInstance().add_motif(base_dir+"/resources/GAAA_tetraloop");

    if(cmd_opts.get_bool("full_path")) {
        auto pb = PathBuilder();
        pb.setup(cmd_opts);
        pb.build();
        
        
        
        return 0;
    }
    
    auto lines = get_lines_from_file(cmd_opts.get_string("path"));
    auto path_points = vectors_from_str(lines[0]);

    auto pf = PathFollower();
    if(cmd_opts.is_filled("mg")) {
        auto lines2 = get_lines_from_file(cmd_opts.get_string("mg"));
        auto mg = std::make_shared<MotifGraph>(lines2[0]);
        
        auto spl = split_str_by_delimiter(lines2[1], " ");
        auto n1 = std::stoi(spl[0]);
        auto end_name = spl[1];
        
        pf.setup(path_points, mg, n1, end_name);
    }
    else {
        pf.setup(path_points);
    }
    
    pf.set_cmd_options(cmd_opts);
    
    if(cmd_opts.get_int("only_one")) {
        auto mt = pf.next();
        
        mt->write_pdbs();
        std::ofstream out;
        out.open("mt_out.top");
        out << mt->topology_to_str();
        out.close();
    }
    
    else {
        
        auto sols = pf.solutions();
        int i = 0;
        for(auto const & s : sols) {
            s->to_pdb("solution."+std::to_string(i)+".pdb", 1);
            i++;
        }
        
        
    }
    
    
}