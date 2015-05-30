//
//  shared_vs_raw.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__shared_vs_raw__
#define __RNAMake__shared_vs_raw__

#include <iostream>
#include <stdio.h>
#include <vector>
#include <memory>

class Node {
public:
    Node(int i) {
        i = i_;
    }
public:
    int i_;
};

class DerivedNode : public Node {
public:
    DerivedNode(int i) : Node(i) {}
};

template <class T>
class Tree {
public:
    Tree() {}
  
    void
    add_data(int i) {
        nodes_1_.push_back(new T(i));
        nodes_2_.push_back(std::shared_ptr<Node>(new Node(i)));
        nodes_3_.push_back(std::unique_ptr<Node>(new Node(i)));

    }
    
    auto
    node(int i ) const {
        return static_cast<T*>(nodes_1_[i]);
    }
    
    std::shared_ptr<T> const &
    node_shared(int i)  {
        node_ = std::static_pointer_cast<T, Node>(nodes_2_[i]);
        return node_;
    }
    
    std::unique_ptr<T> const &
    node_unique(int i)  {
        return std::static_pointer_cast<T, Node>(nodes_2_[i]);
    }
    
private:
    std::vector<Node*> nodes_1_;
    std::vector<std::shared_ptr<Node>> nodes_2_;
    std::vector<std::unique_ptr<Node>> nodes_3_;
    std::shared_ptr<T> node_;
    std::unique_ptr<T> node_2_;
};

#endif /* defined(__RNAMake__shared_vs_raw__) */
