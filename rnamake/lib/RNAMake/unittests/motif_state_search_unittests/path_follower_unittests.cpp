//
//  path_follower_unittests.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/18/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//


#include "resources/resource_manager.h"
#include "motif_data_structures/motif_graph.h"
#include "motif_state_search/path_follower.h"
#include "util/settings.h"
#include "util/file_io.h"
#include "path_follower_unittests.h"


namespace unittests {
namespace motif_state_search {

int
PathFollowerUnittest::test_creation() {
    auto pf = PathFollower();
    return 0;
}
    
int
PathFollowerUnittest::test_build() {
    auto base_dir = String("/Users/josephyesselman/projects/RNAMake/rnamake/lib/RNAMake/apps/mini_ttr");
    ResourceManager::getInstance().add_motif(base_dir+"/resources/GAAA_tetraloop");
    String path = lib_path() + "/unittests/resources/motif_state_search/path_finder.top";
    auto lines = get_lines_from_file(path);
    String path2 = lib_path() + "/unittests/resources/motif_state_search/all_points.str";
    auto lines2 = get_lines_from_file(path2);
    auto path_points = vectors_from_str(lines2[0]);

    auto mg = std::make_shared<MotifGraph>(lines[0]);
    auto spl_1 = split_str_by_delimiter(lines[1], " ");
    auto n1 = std::stoi(spl_1[0]);
    auto name_1 = spl_1[1];
    
    
    auto pf = PathFollower();
    pf.setup(path_points, mg, n1, name_1);
    auto mt = pf.next();
    mg->add_motif_tree(mt, n1, name_1);
    mg->write_pdbs();
    
    
    return 0;
}
    
int
PathFollowerUnittest::run() {
    test_creation();
    test_build();
    return 0;
}
    

}
}