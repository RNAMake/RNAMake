//
//  motif_tree.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/4/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree__
#define __RNAMake__motif_tree__

#include <stdio.h>
#include <map>
#include "types.h"
#include "motif_tree.fwd.h"
#include "motif.h"
#include "uuid.h"

typedef std::map<Uuid, int> UuidIntMap;



class MotifTree {
public:
    MotifTree();
    MotifTree(MotifOP);
    ~MotifTree() {}
    
public:
    
    MotifTreeNodeOP
    add_motif(
        MotifOP const &,
        MotifTreeNodeOP,
        int,
        int);
    
    void
    write_pdbs();

public: //getters
    
    inline
    MotifTreeNodeOPs const &
    nodes() { return nodes_; }
    
private:
    
    MotifTreeNodeOP
    _add_motif(
        MotifTreeNodeOP const &,
        MotifOP const &,
        BasepairOP const &,
        BasepairOP const &,
        Ints const &);

private:
    MotifTreeNodeOPs nodes_;
    MotifTreeNodeOP last_node_;
    
    
};

class MotifTreeNode {
public:
    MotifTreeNode() {}
    MotifTreeNode(
        MotifOP const &,
        int const,
        int const,
        int const);
    
    ~MotifTreeNode() {}

public:
    
    BasepairOPs
    available_ends();
    
    void
    set_end_status(BasepairOP const &, int);
    
    int
    get_end_status(BasepairOP const &);
    
    void
    add_connection(MotifTreeConnection const &);
    
public: //getters:
    
    inline
    MotifOP const &
    motif() const { return motif_; }
    
private:
    int level_, index_, flip_;
    MotifOP motif_;
    UuidIntMap end_status_;
    MotifTreeConnections connections_;
    
};


class MotifTreeConnection {
public:
    MotifTreeConnection() {}
    MotifTreeConnection(
        MotifTreeNodeOP node_1,
        MotifTreeNodeOP node_2,
        BasepairOP end_1,
        BasepairOP end_2,
        int no_overlap=0):
        node_1_( node_1 ),
        node_2_( node_2 ),
        end_1_ ( end_1 ),
        end_2_ ( end_2 ),
        no_overlap_ ( no_overlap )
    {
        node_1_->set_end_status(end_1_, 0);
        node_2_->set_end_status(end_2_, 0);
        node_1_->add_connection(*this);
        node_2_->add_connection(*this);
    }
    
    ~MotifTreeConnection() {}
   
private:
    MotifTreeNodeOP node_1_, node_2_;
    BasepairOP end_1_, end_2_;
    int no_overlap_;

};



#endif /* defined(__RNAMake__motif_tree__) */
