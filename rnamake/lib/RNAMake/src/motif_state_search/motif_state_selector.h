//
//  motif_state_selector.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/11/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_selector__
#define __RNAMake__motif_state_selector__

#include <stdio.h>

//RNAMake Headers
#include "base/types.h"
#include "data_structure/graph/graph.h"
#include "motif/motif_state.h"
#include "motif/motif_state_ensemble.h"
#include "resources/motif_state_sqlite_library.h"
#include "motif_state_search/motif_state_search_node.fwd.h"
#include "motif_state_search/motif_state_search_node.h"


struct MotifStateSelectorNodeData {
    MotifStateSelectorNodeData(
        String const & nname,
        motif::MotifStateOPs const & nmotif_states,
        int nmax_uses=1000,
        int nrequired_uses=0):
    name(nname),
    motif_states(nmotif_states),
    max_uses(nmax_uses),
    required_uses(nrequired_uses)
    {}
    
    String name;
    motif::MotifStateOPs motif_states;
    int max_uses, required_uses;
};

struct MotifStateandType {
    inline
    MotifStateandType():
    motif_state(nullptr),
    type(0)
    {}
    
    inline
    MotifStateandType(
        motif::MotifStateOP const & nmotif_state,
        int ntype):
    motif_state(nmotif_state),
    type(ntype)
    {}
    
    motif::MotifStateOP motif_state;
    int type;
};

class MotifStateandTypes {
public:
    MotifStateandTypes():
    motif_states_and_types_(100),
    pos_(0),
    max_pos_(100)
    {}
    
    ~MotifStateandTypes() {}
    
public: //iterator stuff
    class iterator;
    friend class iterator;
    
    class iterator {
    public:
        inline
        iterator(
            std::vector<MotifStateandType> const & motif_states_and_types,
            int pos,
            int i):
        motif_states_and_types_(motif_states_and_types),
        pos_(pos),
        i_(i)
        {}
        
        iterator operator++() { i_++; return *this; }
        MotifStateandType const & operator*() { return motif_states_and_types_[i_]; }
        bool operator== (iterator const & rhs) const { return i_ == rhs.i_; }
        bool operator!= (iterator const & rhs) const { return i_ != rhs.i_; }
        
    private:
        std::vector<MotifStateandType> motif_states_and_types_;
        int i_;
        int pos_;
    };
    
    iterator begin() { return iterator(motif_states_and_types_, pos_, 0); }
    iterator end()   { return iterator(motif_states_and_types_, pos_, pos_); }

public:
    
    inline
    void
    add(
        motif::MotifStateOP const & ms,
        int type) {
        motif_states_and_types_[pos_].motif_state = ms;
        motif_states_and_types_[pos_].type = type;
        pos_ += 1;
    }
    
    inline
    void
    reset() { pos_ = 0; }
    
    inline
    void
    need_resize(
        int add) {
        if(add + pos_ > max_pos_) {
            max_pos_ = add + pos_ + 1;
            motif_states_and_types_.resize(max_pos_);
        }
    }
    
    inline
    int
    size() { return pos_; }
    
    inline
    std::vector<MotifStateandType> const &
    motif_states_and_types() { return motif_states_and_types_; }
    
    inline
    int
    pos() { return pos_; }
    
private:
    std::vector<MotifStateandType> motif_states_and_types_;
    int max_pos_;
    int pos_;
};

typedef std::shared_ptr<MotifStateSelectorNodeData> MotifStateSelectorNodeDataOP;
//typedef GraphNodeDynamic<MotifStateSelectorNodeDataOP> MotifStateSelectorNodeOP;

class MotifStateSelector {
public:
    MotifStateSelector():
        motif_states_and_types_(MotifStateandTypes()), name_("MotifStateSelector")
    {}
    
    virtual
    ~MotifStateSelector() {}
    
public:
    
    virtual
    void
    add(
        String lib_name = "",
        motif::MotifStateEnsembleOP mse = nullptr,
        motif::MotifOP m = nullptr,
        int max_uses=1000,
        int required_uses=0) {
        
        if(lib_name.length() > 0) {
            auto ms_lib = resources::MotifStateSqliteLibrary(lib_name);
            ms_lib.load_all();
            motif::MotifStateOPs motif_states;
            for(auto const & ms : ms_lib) { motif_states.push_back(ms); }
            auto d = std::make_shared<MotifStateSelectorNodeData>(lib_name, motif_states,
                                                                  max_uses, required_uses);
            graph_.add_data(d);
        }
        
        else if(m != nullptr) {

            auto d = std::make_shared<MotifStateSelectorNodeData>(m->name(),
                                                                  motif::MotifStateOPs { m->get_state() },
                                                                  max_uses,  required_uses);
            graph_.add_data(d);
        }

        else if(mse != nullptr) {
            auto motif_states = motif::MotifStateOPs();
            for(auto const & mem : mse->members()) {
                motif_states.push_back(mem->motif_state);
            }
            auto d = std::make_shared<MotifStateSelectorNodeData>("mse",
                                                                  motif_states,
                                                                  max_uses, required_uses);
            graph_.add_data(d);
        }
    }
    
    int
    size() { return (int)graph_.size(); }

    inline
    String const &
    name() { return name_; }
    
public:
    
    void
    connect(
        String const &,
        String const &);

    void
    connect(
            int,
            int);

    virtual
    MotifStateandTypes const &
    get_children_ms(
            MotifStateSearchNodeOP const & node) {

        motif_states_and_types_.reset();
        //beginning of search start on first node
        if(node->ntype() == -1) {
            motif_states_and_types_.need_resize((int)graph_.get_node(0)->data()->motif_states.size());
            for(auto const & ms : graph_.get_node(0)->data()->motif_states) {
                motif_states_and_types_.add(ms, 0);
            }
        }

        else {
            for(auto const & c : graph_.get_node(node->ntype())->connections()) {
                auto partner =  c->partner(graph_.get_node(node->ntype())->index());
                //did use all that we needed of this type
                //if(partner->data()->max_uses <= node->node_type_usage(partner->index())) {
                //    continue;
                //}
                motif_states_and_types_.need_resize((int)partner->data()->motif_states.size());
                for(auto const & ms : partner->data()->motif_states) {
                    motif_states_and_types_.add(ms, partner->index());
                }
            }
        }

        return motif_states_and_types_;
    }
    
    int
    is_valid_solution(
        MotifStateSearchNodeOP const &);
    
    float
    score(
        MotifStateSearchNodeOP const &);
    
    
protected:
    data_structure::graph::GraphDynamic<MotifStateSelectorNodeDataOP> graph_;
    MotifStateandTypes motif_states_and_types_;
    String name_;
};

class MSS_RoundRobin : public MotifStateSelector {
public:
    MSS_RoundRobin() : MotifStateSelector() {}
    
    ~MSS_RoundRobin() {}
    
public:
    
    void
    add(
        String lib_name = "",
        motif::MotifStateEnsembleOP mse = nullptr,
        motif::MotifOP m = nullptr,
        int max_uses=1000,
        int required_uses=0) {
    
        MotifStateSelector::add(lib_name, mse, m, max_uses, required_uses);
        int i = (int)graph_.size()-1;
        for(auto const & n : graph_) {
            graph_.connect(n->index(), i);
        }
        
    }
    
};

class MSS_HelixFlank : public MotifStateSelector {
public:
    MSS_HelixFlank() : MotifStateSelector() {
        MotifStateSelector::add("ideal_helices_min");
    }
    
    ~MSS_HelixFlank() {}
    
public:
    
    void
    add(
        String lib_name = "",
        motif::MotifStateEnsembleOP mse = nullptr,
        motif::MotifOP m = nullptr,
        int max_uses=1000,
        int required_uses=0) {
        
        MotifStateSelector::add(lib_name, mse, m, max_uses, required_uses);
        int i = (int)graph_.size()-1;
        if(i != 0) {
            graph_.connect(0, i);
        }
    }
    
};

class MSS_Path : public MotifStateSelector {
public:
    MSS_Path() : MotifStateSelector()
    { name_ = "MSS_Path"; }

    ~MSS_Path() {}
public:

    virtual
    MotifStateandTypes const &
    get_children_ms(
            MotifStateSearchNodeOP const & node) {

        motif_states_and_types_.reset();
        //beginning of search start on first node
        if(node->ntype() == -1) {
            motif_states_and_types_.need_resize((int)graph_.get_node(0)->data()->motif_states.size());
            for(auto const & ms : graph_.get_node(0)->data()->motif_states) {
                motif_states_and_types_.add(ms, 0);
            }
        }
        else {
            auto connections =  graph_.get_node(node->ntype())->connections();
            if(connections.size() < 2) { return motif_states_and_types_; }

            if(graph_.size() == node->ntype()+1) {
                return motif_states_and_types_;
            }

            auto partner =  graph_.get_node(node->ntype()+1);
            //std::cout << node->ntype() << " " << partner->index() << std::endl;
            //did use all that we needed of this type
            if(partner->data()->max_uses <= node->node_type_usage(partner->index())) {
                return motif_states_and_types_;
            }
            motif_states_and_types_.need_resize((int)partner->data()->motif_states.size());
            for(auto const & ms : partner->data()->motif_states) {
                motif_states_and_types_.add(ms, partner->index());
            }
        }

        return motif_states_and_types_;
    }

};

typedef std::shared_ptr<MotifStateSelector> MotifStateSelectorOP;

MotifStateSelectorOP
default_selector();



#endif /* defined(__RNAMake__motif_state_selector__) */
































