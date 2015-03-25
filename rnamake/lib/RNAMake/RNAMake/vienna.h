//
//  vienna.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef RNAMake_vienna_h
#define RNAMake_vienna_h

#include "types.h"

typedef struct plist {
    int i;
    int j;
    float p;
    int type;
} plist;

extern "C" float fold(char const*,char*);
extern "C" float cofold(char const*,char*);
extern "C" plist** get_finished_plist(const char *);

struct FoldResult {
public:
    FoldResult() {}
    FoldResult(
        float nfree_energy,
        String nstructure):
    free_energy ( nfree_energy ),
    structure (nstructure) {}
public:
    float free_energy;
    String structure;
};

class Vienna {
public:
    Vienna():
    fold_result_ (FoldResult()) {}
    
    ~Vienna() {}
    
public:
    
    inline
    FoldResult const &
    vfold(
        String const & seq) {
        
        const char * rec_sequence = seq.c_str();
        char * structure = new char[seq.length()];
        float min_en = fold(rec_sequence, structure);
        fold_result_.structure = String(structure);
        fold_result_.free_energy = min_en;
        return fold_result_;
    }

    inline
    FoldResult const &
    vcofold(
        String const & seq) {
        
        const char * rec_sequence = seq.c_str();
        char * structure = new char[seq.length()];
        float min_en = cofold(rec_sequence, structure);
        fold_result_.structure = String(structure);
        fold_result_.free_energy = min_en;
        return fold_result_;
    }
    
    inline
    plist** const &
    dotplot(
        String const & seq) {
        
        //std::cout << pl1[4]->p << std::endl;
        //std::cout << probs[0] << std::endl;
        if(used_ == 1) { delete pl1_; }
        
        pl1_ = get_finished_plist(seq.c_str());
        used_ = 1;
        //std::cout << "start" << std::endl;
        
        for(int i = 0; i < seq.length()*seq.length(); i++) {
            //if(pl1_[i]->p < 0.0001) { continue; }
            //std::cout << pl1_[i]->i << " " << pl1_[i]->j << " " << pl1_[i]->p << std::endl;
        }
        
        //delete pl1;
        //std::cout << test << std::endl;
        return pl1_;
    }

    
private:
    FoldResult fold_result_;
    plist** pl1_;
    int used_;
    
};

#endif
