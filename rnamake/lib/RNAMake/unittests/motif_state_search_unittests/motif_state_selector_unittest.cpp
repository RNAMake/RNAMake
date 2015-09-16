//
//  motif_state_selector_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/13/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_state_selector_unittest.h"
#include "resources/resource_manager.h"
#include "motif_state_search/motif_state_selector.h"
#include "motif_state_search/motif_state_search_node.h"

int
MotifStateSelectorUnittest::test_creation() {
    MotifStateSelector mss;
    mss.add("twoway");
    return 1;
}

int
MotifStateSelectorUnittest::test_get_children_ms() {
    MotifStateSelector mss;
    mss.add("unique_twoway");
    auto ms = ResourceManager::getInstance().get_state("HELIX.IDEAL.2");
    auto node = std::make_shared<MotifStateSearchNode>(ms, nullptr, -1, -1);
    auto children = mss.get_children_ms(node);
    return 1;
}

int
MotifStateSelectorUnittest::test_default_selector() {
    auto selector = default_selector();
    auto ms = ResourceManager::getInstance().get_state("HELIX.IDEAL.2");
    auto node = std::make_shared<MotifStateSearchNode>(ms, nullptr, -1, -1);
    auto children = selector->get_children_ms(node);
    //std::cout << children.size() << std::endl;
    return 1;
}


int
MotifStateSelectorUnittest::run() {
    if(test_creation() == 0)         { std::cout << "test_creation failed" << std::endl; }
    if(test_get_children_ms() == 0)  { std::cout << "test_get_children_ms failed" << std::endl; }
    if(test_default_selector() == 0) { std::cout << "test_default_selector failed" << std::endl; }
    return 1;
}

