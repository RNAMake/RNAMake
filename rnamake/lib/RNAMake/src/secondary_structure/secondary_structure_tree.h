//
//  secondary_structure_tree.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__secondary_structure_tree__
#define __RNAMake__secondary_structure_tree__

#include <stdio.h>

#include "data_structure/tree/tree.h"
#include "secondary_structure/motif.h"
#include "secondary_structure/pose.h"

namespace sstruct {

class SecondaryStructureTree {
public:
    SecondaryStructureTree():
    tree_(TreeStatic<MotifOP>())
    {}
    
    
    ~SecondaryStructureTree() {}
    
public: //iterators
    
    typedef typename TreeStatic<MotifOP>::iterator iterator;
    typedef typename TreeStatic<MotifOP>::const_iterator const_iterator;
    
    iterator begin() { return tree_.begin(); }
    iterator end()   { return tree_.end(); }
    
    const_iterator begin() const { return tree_.begin(); }
    const_iterator end()   const { return tree_.end(); }
    
public:
    size_t
    size() { return tree_.size(); }
    
    int
    add_motif(
        MotifOP const & m,
        int parent_index = -1,
        int parent_end_index = -1);
    
private:
    TreeStatic<MotifOP> tree_;
};

typedef std::shared_ptr<SecondaryStructureTree> SecondaryStructureTreeOP;

SecondaryStructureTreeOP
tree_from_pose(PoseOP const &);
    

    
}

#endif /* defined(__RNAMake__secondary_structure_tree__) */
