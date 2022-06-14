#ifndef __tools_builder_graph
#define __tools_builder_graph


#include "data_structure/graph/graph.h"
#include "util/motif_type.h"

class BuilderGraph {
public:
    
    BuilderGraph() {
        g_ = data_structure::graph::GraphStatic<util::MotifType>();
    }
    
    size_t
    size() { return g_.size(); }
    
    int
    add_node(
        util::MotifType type,
        int ends,
        int parent_index = -1,
        int parent_end_index = -1) {
        
        auto parent = g_.last_node();
        
        if(parent_index != -1) {
            parent = g_.get_node(parent_index);
        }
        
        if(parent == nullptr) {
            return g_.add_data(type, -1, -1, -1, ends);
        }
        
        Indexes avail_pos;
        if(parent_end_index == -1) {
            avail_pos = parent->available_children_pos();
        }
        else {
            avail_pos = Indexes{parent_end_index};
        }
        
        if(avail_pos.size() == 0) { return -1; }
        
        for(auto const & p : avail_pos) {
            if(p == 0) { continue; }
            
            return g_.add_data(type, parent->index(), p, 0, ends);
        }
        
        return -1;
    }
    
    inline
    data_structure::graph::GraphNodeOP<util::MotifType> const &
    get_node(int i) { return g_.get_node(i); }
    
public: //iterators
    
    typedef typename data_structure::graph::GraphStatic<util::MotifType>::iterator iterator;
    typedef typename data_structure::graph::GraphStatic<util::MotifType>::const_iterator const_iterator;
    
    iterator begin() { return g_.begin(); }
    iterator end()   { return g_.end(); }
    
    const_iterator begin() const { return g_.begin(); }
    const_iterator end()   const { return g_.end(); }
    
    
private:
    data_structure::graph::GraphStatic<util::MotifType> g_;
};


typedef std::shared_ptr<BuilderGraph> BuilderGraphOP;

BuilderGraphOP
helix_and_two_way() {
    auto g = std::make_shared<BuilderGraph>();
    g->add_node(util::MotifType::HELIX, 2);
    g->add_node(util::MotifType::TWOWAY, 2);
    
    return g;
}


BuilderGraphOP
helix_and_two_way_and_hairpin() {
    auto g = std::make_shared<BuilderGraph>();
    g->add_node(util::MotifType::HELIX, 2);
    g->add_node(util::MotifType::TWOWAY, 2);
    g->add_node(util::MotifType::HELIX, 2);
    g->add_node(util::MotifType::HAIRPIN, 1);
    return g;

}




#endif
