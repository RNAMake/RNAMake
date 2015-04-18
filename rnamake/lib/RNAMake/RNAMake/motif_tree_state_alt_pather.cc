//
//  motif_tree_state_alt_pather.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_state_alt_pather.h"
#include "settings.h"
#include "FileIO.h"
#include "cartesian_product.h"

MotifTreeStateAltPather::MotifTreeStateAltPather() {
    
    String lib_path = resources_path() + "/precomputed/motif_tree_states/TWOWAY_all.new.me";
    mts_lib_ = MotifTreeStateLibrary(lib_path);
    all_trees_ = std::vector<MotifTreeStateTree>();
    sim_dict_ = std::map<String, Strings>();
    String path = resources_path() + "sim_list";
    String line;
    std::ifstream in;
    in.open(path);
    
    int i = -1;
    while (in.good()) {
        getline(in, line);
        if(line.length() < 5) { break; }
        i++;
        if(i == 0) { continue;}
        
        Strings spl = split_str_by_delimiter(line, " ");
        String key = spl[1];
        Strings alts (spl.size()-4);
        
        for(int i = 4; i < spl.size(); i++) {
            alts[i-4] = spl[i];
        }
        
        sim_dict_[key] = alts;
    }
}

std::vector<MotifTreeStateTree> const &
MotifTreeStateAltPather::get_alt_paths(
    MotifTreeStateTree const & mtst) {
    
    std::vector<MotifTreeStateOPs> mts_alts;
    MotifTreeStateOPs all_mts;
    MotifTreeStateOP mts;
    all_trees_.resize(0);
    
    Strings name_spl, alt_motifs;
    String rest, mname, mts_name;
    int i = -1;
    for(auto const & n : mtst.nodes()) {
        i++;
        if(i == 0) { continue; }
        name_spl = split_str_by_delimiter(n->mts()->name(), "-");
        mname = name_spl[0];
        rest = "";
        for(int j = 1; j < name_spl.size(); j++) {
            rest += name_spl[j];
            if(j < name_spl.size()-1 ) {
                rest += "-";
            }
        }
        
        if (sim_dict_.find(mname) != sim_dict_.end()) {
            alt_motifs = sim_dict_[mname];
            MotifTreeStateOPs all_mts;
            for(auto const & alt_name : alt_motifs) {
                mts_name = alt_name + "-" + rest;
                mts = mts_lib_.get_state_no_error(mts_name);
                if(mts == NULL) { continue; }
                all_mts.push_back(mts);
            }
            mts_alts.push_back(all_mts);
        }
        
        else {
            MotifTreeStateOPs all_mts(1);
            all_mts[0] = n->mts();
            mts_alts.push_back(all_mts);
        }
        
    }
   
    CartesianProduct<MotifTreeStateOP> product(mts_alts);
    MotifTreeStateOPs c;
    MotifTreeStateNodeOP parent, node;
    int end_index, fail, count = 0;
    
    
    while(!product.end()) {
        MotifTreeStateTree mtst2;
        c = product.next();
        i = 1;
        fail = 0;
        for(auto const & mts : c) {
            parent = mtst2.nodes()[ mtst.nodes()[i]->parent()->index() ];
            end_index = mtst.nodes()[i]->parent_end_index();
            
            node = mtst2.add_state(mts, parent, end_index);
            if(node == NULL) { fail = 1; break; }
            i++;
        }
        
        if(fail) { continue;}
        
        all_trees_.push_back(mtst2);

    }
    
    return all_trees_;
        
    
    
    
    
}