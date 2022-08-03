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

namespace secondary_structure {

class SecondaryStructureTree {
public:
  SecondaryStructureTree()
      : tree_(data_structure::tree::TreeStatic<MotifOP>()) {}

  ~SecondaryStructureTree() {}

public: // iterators
  typedef typename data_structure::tree::TreeStatic<MotifOP>::iterator iterator;
  typedef typename data_structure::tree::TreeStatic<MotifOP>::const_iterator
      const_iterator;

  iterator begin() { return _tree.begin(); }
  iterator end() { return _tree.end(); }

  const_iterator begin() const { return _tree.begin(); }
  const_iterator end() const { return _tree.end(); }

public:
  size_t size() { return _tree.size(); }

  int add_motif(MotifOP const &m, int parent_index = -1,
                int parent_end_index = -1);

private:
  data_structure::tree::TreeStatic<MotifOP> _tree;
};

typedef std::shared_ptr<SecondaryStructureTree> SecondaryStructureTreeOP;

SecondaryStructureTreeOP tree_from_pose(PoseOP const &);

} // namespace secondary_structure

#endif /* defined(__RNAMake__secondary_structure_tree__) */
