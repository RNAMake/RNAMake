#ifndef __tools_build_motif_graph
#define __tools_build_motif_graph

#include <map>

#include "builder_graph.hpp"

#include "motif_data_structure/motif_graph.h"
#include "resources/motif_sqlite_library.h"


class MotifGraphBuilder {
public:
    MotifGraphBuilder() {
        g_ = helix_and_two_way();
        mlibs_ = std::map<String, resources::MotifSqliteLibraryOP>();
    }
    
    MotifGraphBuilder(
        BuilderGraphOP const & g) {
        
        g_ = g;
        mlibs_ = std::map<String, resources::MotifSqliteLibraryOP>();
    }
    
    
public:
    motif_data_structure::MotifGraphOP
    build(int repeat=1) {
        auto mg = std::make_shared<motif_data_structure::MotifGraph>();
        
        for(int i = 0; i < repeat; i++) {
            auto index_map = std::map<int, int>();
            int j = 0, pos = 0, parent_index = -1, parent_end_index = -1;
            for(auto const & n : *g_) {
                auto lib_name = _get_lib_name(n->data());
                auto mlib = _get_lib(lib_name);
                
                if(j != 0) {
                    parent_index = n->parent()->index();
                    parent_end_index = n->parent_end_index();
                    
                    parent_index = index_map[parent_index];
                }
                
                pos = _add_motif_to_graph(mg, mlib, parent_index, parent_end_index);
                
                if(pos == -1) { return mg; }
                index_map[n->index()] = pos;
                j++;
                
            }
        }
        
        return mg;
    }
    
private:
    String
    _get_lib_name(util::MotifType const & type) {
        if     (type == util::MotifType::HELIX)   { return "ideal_helices"; }
        else if(type == util::MotifType::TWOWAY)  { return "twoway"; }
        else if(type == util::MotifType::NWAY)    { return "nway"; }
        else if(type == util::MotifType::HAIRPIN) { return "hairpin"; }
        else if(type == util::MotifType::TCONTACT){ return "tcontact"; }
        else { throw std::runtime_error("cannot convert type"); }
        
    }
    
    resources::MotifSqliteLibraryOP const &
    _get_lib(String const & name) {
        if(mlibs_.find(name) == mlibs_.end()) {
            mlibs_[name] = std::make_shared<resources::MotifSqliteLibrary>(name);
        }
        
        return mlibs_[name];
    }
    
    inline
    int
    _add_motif_to_graph(
        motif_data_structure::MotifGraphOP & mg,
        resources::MotifSqliteLibraryOP const & mlib,
        int parent_index,
        int parent_end_index) {
        
        int pos = -1, count = 0;
        while(pos == -1) {
            auto m = mlib->get_random();
            
            pos = mg->add_motif(m, parent_index, parent_end_index);
            
            count++;
            if(count > 100) {
                return 0;
            }
        }
        
        return pos;
    }
    
    
private:
    BuilderGraphOP g_;
    std::map<String, resources::MotifSqliteLibraryOP> mlibs_;
    
};


#endif














