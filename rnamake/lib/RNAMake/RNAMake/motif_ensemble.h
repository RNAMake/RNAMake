//
//  motif_ensemble.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/6/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_ensemble__
#define __RNAMake__motif_ensemble__

#include <map>
#include <fstream>
#include <random>
#include <stdio.h>
#include <vector>
#include "types.h"
#include "settings.h"
#include "FileIO.h"
#include "motif_tree_state.h"
#include "motif_tree_state_library.h"


struct MotifState {
public:
    MotifState(
        MotifTreeStateOP const & nmts,
        float const & npopulation):
        mts(nmts),
        population(npopulation){
        }
    
    ~MotifState() {}
    
public:
    MotifTreeStateOP mts;
    float population;

};

struct less_than_key
{
    inline bool operator() (MotifState const & struct1, MotifState const & struct2)
    {
        return (struct1.population < struct2.population);
    }
};

typedef std::vector<MotifState> MotifStates;

class MotifEnsemble {
public:
    MotifEnsemble(
        String const lib_path = "",
        int const start_index = -1,
        int const flip_direction = -1
    ) {
        
        motif_states_ = MotifStates();
        if(lib_path.length() < 2) { return; }
        
        String prediction_path = resources_path()+"/prediction/";
        String mts_lib_path = prediction_path+lib_path+".new.me";
        String pop_path = prediction_path+lib_path+".pop";
        if (!file_exists(mts_lib_path)) {
            mts_lib_path = lib_path+".new.me";
            pop_path = lib_path+".pop";
        }
        String line;
        StringFloatMap pop;
        
        std::ifstream input;
        input.open(pop_path.c_str());
        while ( input.good() ) {
            getline(input, line);
            if ( line.length() < 10) { continue; }
            Strings spl = split_str_by_delimiter(line, " ");
            pop[ spl[0] ] = std::stof(spl.back());
        }
        
        input.close();
        MotifTreeStateLibrary mts_lib ( mts_lib_path );
        for (auto const & mts : mts_lib.motif_tree_states()) {
            if (start_index != -1 && start_index != mts->start_index()) { continue; }
            if (flip_direction != -1 && flip_direction != mts->flip()) { continue; }
            NameElements name_elements = parse_db_name(mts->name());
            float m_pop = pop [ name_elements.motif_name ];
            MotifState motif_state ( mts, m_pop);
            motif_states_.push_back(motif_state);
        }
        std::sort(motif_states_.begin(),motif_states_.end(),less_than_key());
        
    }
    
    ~MotifEnsemble() {}

public:
    void
    add_motif_state(MotifState const & ms) {
        motif_states_.push_back(ms);
    }
    
    MotifState const &
    get_state(String const &) const;
    
    MotifState const &
    get_random_state() const;
    
    MotifState const &
    get_state(int const & pos) const { return motif_states_[pos]; }
    
public: // getters
    
    inline
    MotifStates const &
    motif_states() const { return motif_states_; }
    
private:
    MotifStates motif_states_;

    
    
};


#endif /* defined(__RNAMake__motif_ensemble__) */
