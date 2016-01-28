//
//  follow_path.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/21/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "follow_path.h"
#include "util/file_io.h"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_graph.h"
#include "motif_state_search/motif_state_search.h"
#include "motif_state_search/path_follower.h"


CommandLineOptions
parse_command_line(
    int argc,
    const char ** argv) {
    
    CommandLineOptions cl_opts;
    auto search = MotifStateSearch();
    cl_opts.add_options(search.options());
    cl_opts.add_option("path", String(""), OptionType::STRING, true);
    cl_opts.add_option("mg", String(""), OptionType::STRING, true);
    cl_opts.parse_command_line(argc, argv);

    return cl_opts;
    
}


int main(int argc, const char * argv[]) {
    auto cmd_opts = parse_command_line(argc, argv);
    
    //load TTR
    auto base_dir = String("/Users/josephyesselman/projects/RNAMake/rnamake/lib/RNAMake/apps/mini_ttr");
    ResourceManager::getInstance().add_motif(base_dir+"/resources/GAAA_tetraloop");
    
    auto lines = get_lines_from_file(cmd_opts.get_string("path"));
    auto lines2 = get_lines_from_file(cmd_opts.get_string("mg"));
    auto path_points = vectors_from_str(lines[0]);
    auto mg = std::make_shared<MotifGraph>(lines2[0]);
    
    auto spl = split_str_by_delimiter(lines2[1], " ");
    auto n1 = std::stoi(spl[0]);
    auto end_name = spl[1];
    
    auto pf = PathFollower();
    pf.setup(path_points, mg, n1, end_name);
    auto mt = pf.next();
    
    mt->write_pdbs();
    std::ofstream out;
    out.open("mt_out.top");
    out << mt->topology_to_str();
    out.close();
    
}