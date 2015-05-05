//
//  motif_tree_node.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_node__
#define __RNAMake__motif_tree_node__

#include <stdio.h>

//RNAMake Headers
#include "structure/basepair.h"
#include "motif/motif_tree.fwd.h"
#include "motif/motif_tree_node.h"
#include "motif/motif.h"

typedef std::map<String, int> UuidIntMap;

class MotifTreeNode {
public:
    MotifTreeNode() {}
    MotifTreeNode(
        MotifOP const &,
        int const,
        int const,
        int const);
    
    ~MotifTreeNode();
public:
    
    BasepairOPs
    available_ends();
    
    void
    set_end_status(BasepairOP const &, int);
    
    int
    get_end_status(BasepairOP const &);
    
    void
    add_connection(MotifTreeConnectionOP const &);
    
public: //getters:
    
    inline
    MotifOP const &
    motif() const { return motif_; }
    
    inline
    int const &
    level() const { return level_; }
    
    inline
    int const &
    index() const { return index_; }
    
    inline
    int const &
    flip() const { return flip_; }
    
    inline
    MotifTreeConnectionOPs const &
    connections() { return connections_; }
    
private:
    int level_, index_, flip_;
    MotifOP motif_;
    UuidIntMap end_status_;
    MotifTreeConnectionOPs connections_;
    
};


class MotifTreeConnection {
public:
    MotifTreeConnection() {}
    MotifTreeConnection(
        MotifTreeNodeOP const & node_1,
        MotifTreeNodeOP const & node_2,
        BasepairOP const & end_1,
        BasepairOP const & end_2,
        int no_overlap=0):
    node_1_( node_1 ),
    node_2_( node_2 ),
    end_1_ ( end_1 ),
    end_2_ ( end_2 ),
    no_overlap_ ( no_overlap )
    {
        node_1_->set_end_status(end_1_, 0);
        node_2_->set_end_status(end_2_, 0);
    }
    
    ~MotifTreeConnection() {
        node_1_ = nullptr;
        node_2_ = nullptr;
        end_1_  = nullptr;
        end_2_  = nullptr;
    
    }
    
public:
    
    bool
    operator == (MotifTreeConnection const & other) const  {
        return (node_1_->index() == other.node_1_->index() && node_2_->index() == other.node_2_->index());
    }
    
    bool
    operator < (MotifTreeConnection const & other) const  {
        return end_1_->uuid().s_uuid().compare(other.end_1_->uuid().s_uuid());
    }
    
public:
    inline
    void
    disconnect() {
        if(node_1_.get() != NULL) {
            node_1_->set_end_status(end_1_, 1);
            node_2_->set_end_status(end_2_, 1);
        }
        node_1_ = nullptr;
        node_2_ = nullptr;
        end_1_ = nullptr;
        end_2_ = nullptr;

    }
    
    inline
     MotifTreeNodeOP const &
     partner(MotifTreeNodeOP const & node) const {
         if     (node == node_1_)  { return node_2_; }
         else if(node == node_2_)  { return node_1_; }
         else { throw "node is not in connection object cannot call partner"; }
     }
     
     inline
     BasepairOP const &
     motif_end(MotifTreeNodeOP const & node) const {
         if     (node == node_1_)  { return end_1_; }
         else if(node == node_2_)  { return end_2_; }
         else { throw "node is not in connection object cannot call motif_end"; }
     
     }
    
private:
    MotifTreeNodeOP node_1_, node_2_;
    BasepairOP end_1_, end_2_;
    int no_overlap_;
    
};


#endif /* defined(__RNAMake__motif_tree_node__) */
