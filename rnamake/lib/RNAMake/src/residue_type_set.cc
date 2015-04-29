//
//  residue_type_set.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <dirent.h>
#include <fstream>
#include "residue_type_set.h"
#include "settings.h"
#include "types.h"

ResidueTypeSet::ResidueTypeSet() {
    residue_types_ = ResidueTypes();
    String path = resources_path() + "/residue_types/";
    DIR *pDIR;
    struct dirent *entry;
    pDIR=opendir(path.c_str());
    while ((entry = readdir(pDIR)) != NULL) {
        String fname ( entry->d_name );
        if(fname.length() < 4) { continue; }
        
        String name = _get_rtype_name(path + fname);
        StringIntMap atom_map = _get_atom_map_from_file(path + fname);
        ResidueType rtype ( name, atom_map);
        residue_types_.push_back(rtype);
    }
    delete pDIR;
    delete entry;
    
}

String
ResidueTypeSet::_get_rtype_name(
    String const & fname) {
    Strings name_spl = split_str_by_delimiter(fname, "/");
    String type_file_name = name_spl.back();
    std::size_t pos = type_file_name.find(".");
    return type_file_name.substr(0,pos);
    
}

StringIntMap
ResidueTypeSet::_get_atom_map_from_file(
    String const & fname) {
    
    std::string line;
    std::ifstream input;
    input.open(fname);
    getline(input, line);
    input.close();
    Strings atom_names = split_str_by_delimiter(line, " ");
    int i = 0;
    StringIntMap atom_map;
    for(auto const & name : atom_names) {
        atom_map[name] = i;
        i++;
    }
    return atom_map;
}

ResidueType
const &
ResidueTypeSet::get_rtype_by_resname(
    String const & resname) const {
    
    for (auto const & restype : residue_types_) {
        if(restype.match_name(resname)) { return restype; }
    }
    
    throw("(cannot find residue with name :"+resname);
}



