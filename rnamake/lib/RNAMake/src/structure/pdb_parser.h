//
//  pdb_parser.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/12/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__pdb_parser__
#define __RNAMake__pdb_parser__

#include <stdio.h>

//RNAMake Headers
#include "structure/residue.h"
#include "structure/residue_type_set.h"

class PDBParser {
public:
    PDBParser():
    residues_(ResidueOPs()),
    rts_(ResidueTypeSet())
    {}
    
    ~PDBParser() {}
    
public:
    
    ResidueOPs const &
    parse(String const & pdb_file,
          int protein=0,
          int rna=1);
    

private:
    ResidueOPs residues_;
    ResidueTypeSet rts_;
    
};

#endif /* defined(__RNAMake__pdb_parser__) */
