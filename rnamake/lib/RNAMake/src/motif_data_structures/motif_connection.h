//
//  motif_connection.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/9/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_connection__
#define __RNAMake__motif_connection__

#include <stdio.h>

#include "base/types.h"

class MotifConnection {
public:
    MotifConnection() {}
    
    inline
    MotifConnection(
        int i,
        int j,
        String const & name_i,
        String const & name_j):
    i_(i),
    j_(j),
    name_i_(name_i),
    name_j_(name_j)
    {}
    
    inline
    MotifConnection(
        MotifConnection const & mc):
    i_(mc.i_),
    j_(mc.j_),
    name_i_(mc.name_i_),
    name_j_(mc.name_j_)
    {}
    
    ~MotifConnection() {}
    
public: // getters
    inline
    int
    i() { return i_; }
    
    inline
    int
    j() { return j_; }
    
    inline
    String const &
    name_i() { return name_i_; }
    
    inline
    String const &
    name_j() { return name_j_; }
    
public:
    
    String
    to_str();

    
private:
    int i_, j_;
    String name_i_, name_j_;
    
    
};

typedef std::shared_ptr<MotifConnection> MotifConnectionOP;
typedef std::vector<MotifConnectionOP>   MotifConnectionOPs;

class MotifConnections {
public:
    MotifConnections();
    
    MotifConnections(
        MotifConnections const &);
    
    ~MotifConnections() {}
    
public: //iterators
    
    typedef typename MotifConnectionOPs::iterator iterator;
    typedef typename MotifConnectionOPs::const_iterator const_iterator;
    
    iterator begin() { return connections_.begin(); }
    iterator end()   { return connections_.end(); }
    
    const_iterator begin() const { return connections_.begin(); }
    const_iterator end()   const { return connections_.end(); }
    
public:
    
    inline
    size_t
    size() const { return connections_.size(); }
    
    void
    add_connection(
        int i,
        int j,
        String const & name_i,
        String const & name_j) {
        
        connections_.push_back(std::make_shared<MotifConnection>(i, j, name_i, name_j));
    }
    
    void
    remove_connections_to(
        int index) {
        
        int pos = 0;
        while(pos < connections_.size()) {
            if(index == connections_[pos]->i()) {
                connections_.erase(connections_.begin()+pos);
                pos--;
            }
            else if(index == connections_[pos]->j()) {
                connections_.erase(connections_.begin()+pos);
                pos--;
            }
            pos++;
            
        }
        
    }
    
    bool
    in_connection(
        int index,
        String const & name) const {
        
        for(auto const & c : connections_) {
            if(index == c->i() && name == c->name_i()) { return true; }
            if(index == c->j() && name == c->name_j()) { return true; }
        }
        return false;
        
    }
    

private:
    MotifConnectionOPs connections_;
};


#endif /* defined(__RNAMake__motif_connection__) */




















