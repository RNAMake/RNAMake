//
//  residue_type.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "residue_type.h"
#include "types.h"
#include "atom.h"

ResidueType::ResidueType(
    String const & name,
    StringIntMap const & atom_map):
    name_ (name),
    atom_map_ ( atom_map) {
    
    alt_names_ = Strings();
    alt_names_.push_back(name_.substr(0,1));
    alt_names_.push_back("r" + name_.substr(0,1));
    alt_names_.push_back("D" + name_.substr(0,1));
    extend_res_specific_altnames();
        
    atom_alt_names_ = StringStringMap();
    atom_alt_names_["O1P"] = "OP1";
    atom_alt_names_["O2P"] = "OP2";
    
}

String
ResidueType::get_correct_atom_name(
    Atom const & a) const {

    StringStringMap::const_iterator iter(atom_alt_names_.find(a.name()));
    if(iter != atom_alt_names_.end()) {
        return iter->second;
    }
    else {
        return "";
    }
}

int
ResidueType::match_name(
    String const & name) const{
    if (name_.compare(name) == 0 ) { return 1; }
    for (auto const & s : alt_names_) {
        if (s.compare(name) == 0) { return 1; }
    }
    return 0;
}

void
ResidueType::extend_res_specific_altnames() {
    // There has to be a better way to do this
    Strings alt_names;
    if      ( name_.compare("GUA") == 0 ) {
        alt_names = split_str_by_delimiter("MIA GDP GTP M2G 1MG 7MG G7M QUO I YG", " ");
    }
    else if ( name_.compare("ADE") == 0 ) {
        alt_names = split_str_by_delimiter("A23 3DA 1MA 12A AET 2MA", " ");
    }
    else if ( name_.compare("URA") == 0 ) {
        alt_names = split_str_by_delimiter("PSU H2U 5MU 4SU 5BU 5MC U3H 2MU 70U BRU DT", " ");
    }
    else {
        alt_names = split_str_by_delimiter("CBR CCC", " ");
    }
    
    for ( auto const & name : alt_names_ ){
        alt_names_.push_back(name);
    }
    
}






