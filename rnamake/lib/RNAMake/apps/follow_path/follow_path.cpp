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
#include "util/settings.h"

void
PathBuilder::setup_options() {
    options_.add_option("path", String(""), OptionType::STRING);
    options_.add_option("mg", String(""), OptionType::STRING);
    options_.add_option("solutions", 1, OptionType::INT);
    options_.add_option("verbose", true, OptionType::BOOL);
    options_.lock_option_adding();
    update_var_options();
}

void
PathBuilder::setup(
    CommandLineOptions const & cmd_opts) {
    
    for(auto const & opt: cmd_opts) {
        if(! has_option(opt->name())) { continue; }
        if(! opt->filled()) { continue; }
        
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
    auto end_bp = BasepairStateOP(nullptr);
    if(get_string_option("mg").size() > 0) {
        auto lines2 = get_lines_from_file(get_string_option("mg"));
        mg = std::make_shared<MotifGraph>(lines2[0]);
        
        auto spl = split_str_by_delimiter(lines2[1], " ");
        n1 = std::stoi(spl[0]);
        end_name = spl[1];
        
        if(lines2.size() > 3) {
            spl = split_str_by_delimiter(lines2[2], " ");
            end_bp = mg->get_end(std::stoi(spl[0]), spl[1])->state();
        }
    }
    
    auto segments = _get_segments(path_points);
    auto pathes = _get_sub_pathes(segments);
    
    auto mt = graph_to_tree(mg, mg->get_node(0));
    auto mst = std::make_shared<MotifStateTree>(mt);
    
    auto pf = PathFollower();
    pf.setup(pathes[0], mg, n1, end_name);
    std::cout << options_.get_int("solutions") << std::endl;
    pf.set_option_value("max_pathes", options_.get_int("solutions"));
    auto sol_org = pf.solutions();
    
    nodes_ = PathBuilderNodes();
    for(int i = 0; i < sol_org.size(); i++) {
        auto sol_mst = sol_org[i]->to_mst();
        auto n = PathBuilderNode { std::make_shared<MotifStateTree>(*mst), 0, 0, 0 };
        n.mst->add_mst(sol_mst, n1, -1, end_name);
        n.add_score(sol_mst, sol_org[i]->score());
        nodes_.push_back(n);
    }

    std::sort(nodes_.begin(), nodes_.end(), PathBuilderNode_LessThanKey());
    
    if(options_.get_bool("verbose")) {
        std::cout << "ROUND 0: SIZE=" << nodes_.size() << " MIN SCORE=" << nodes_[0].path_score;
        std::cout << std::endl;
    }
    

    //auto path = pathes.back();
    //pathes.pop_back();
    //for(auto const & p : path) {
    //    pathes.back().push_back(p);
    //}
    
    auto new_nodes = PathBuilderNodes();
    for(int level = 1; level < pathes.size(); level++) {
        new_nodes = PathBuilderNodes();
        
        for(auto const & n : nodes_) {
            
            pf = PathFollower();
            if(level == pathes.size()-1 && end_bp != nullptr) {
                pf.setup(pathes[level], n.mst->last_node()->data()->cur_state->end_states()[1],
                         end_bp, n.mst->centers());
            }

            else {
                pf.setup(pathes[level], n.mst->last_node()->data()->cur_state->end_states()[1],
                         n.mst->centers());
            }
            int num = options_.get_int("solutions")/10;
            if(num > 1) { num = 1; }
            pf.set_option_value("max_pathes", num);
            
            
            auto sols = pf.solutions();
            
            for(auto const & sol : sols) {
                auto new_n = PathBuilderNode(n);
                new_n.add_solution(sol);
                new_nodes.push_back(new_n);
            }
            
        }
        
        std::sort(new_nodes.begin(), new_nodes.end(), PathBuilderNode_LessThanKey());
        if(new_nodes.size() > options_.get_int("solutions")) {
            nodes_ = PathBuilderNodes(&new_nodes[0], &new_nodes[options_.get_int("solutions")]);
        }
        else {
            nodes_ = new_nodes;

        }
        
        if(nodes_.size() == 0) {
            std::cout << "full run, ran out of options" << std::endl;
            return;
        }
            
        if(options_.get_bool("verbose")) {
            std::cout << "ROUND " << level << ": SIZE=" << nodes_.size();
            std::cout <<" MIN SCORE=" << nodes_[0].path_score << std::endl;
        }
        
        
    }
    
    std::ofstream out2;
    out2.open("solutions.top");
    
    int i = 0;
    for(auto const & n : nodes_) {
        n.mst->to_motif_tree()->to_pdb("solution."+std::to_string(i)+".pdb", 1);
        out2 << n.mst->topology_to_str() << std::endl;
        i += 1;
    }
    
    out2.close();
    
    
    auto final_mst = std::make_shared<MotifStateTree>();
    for(auto const & n : *nodes_[0].mst) {
        if(n->level() == 0) { continue; }
        final_mst->add_state(n->data()->cur_state, -1, -1, "", true);
    }

    std::ofstream out;
    out.open("mt_out.top");
    out << final_mst->topology_to_str();
    out.close();
    
 
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
            for(int j = 0; j < half-1; j++) { p.push_back(segments[i][j]); }
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
    auto ttr_dir = String(base_dir()+"/rnamake/lib/RNAMake/apps/mini_ttr");
    ResourceManager::getInstance().add_motif(ttr_dir+"/resources/GAAA_tetraloop");

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
