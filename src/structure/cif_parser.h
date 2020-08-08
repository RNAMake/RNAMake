//
//  cif_parser.h
//  RNAMake
//
//  Created by Chris Jurich on Aug 2 2020
//  Copyright (c) 2020 Chris Jurich. All rights reserved.
//

#ifndef __RNAMake__cif_parser__
#define __RNAMake__cif_parser__

#include <stdio.h>
#include <exception>
//RNAMake Headers
#include "structure/residue.h"
#include "structure/residue_type_set.h"


namespace structure {

class CIFParseError : public std::runtime_error {
    public:
        CIFParseError(
                String& msg
                ) : std::runtime_error(msg) {}
};


class CIFParser {
public:
    CIFParser() :
            residues_(ResidueOPs()),
            rts_(ResidueTypeSet()) {}

    ~CIFParser() {}

public:

    ResidueOPs const &
    parse(
            String const & pdb_file,
            int protein = 0,
            int rna = 1,
            int others = 0);

private:

    ResidueType
    _get_new_residue_type(
            String const &,
            Strings const &);

private:
    ResidueOPs residues_;
    ResidueTypeSet rts_;

private:
    std::unordered_set<int> 
    _get_entity_ids(
                String const &
                    );

};

}

#endif /* defined(__RNAMake__cif_parser__) */
