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

extern "C" float fold(char const*,char*);

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

private:
    FoldResult fold_result_;
    
};

#endif
