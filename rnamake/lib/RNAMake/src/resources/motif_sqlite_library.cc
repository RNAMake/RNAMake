//
//  motif_sqlite_library.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "resources/motif_sqlite_library.h"
#include "structure/residue_type_set_manager.h"


StringStringMap
MotifSqliteLibrary::get_libnames() {
    StringStringMap libnames;
    
    libnames["ideal_helices"]  = "/motif_libraries_new/ideal_helices.db";
    libnames["twoway"]         = "/motif_libraries_new/twoway.db";
    libnames["tcontact"]       = "/motif_libraries_new/tcontact.db";
    libnames["hairpin"]        = "/motif_libraries_new/hairpin.db";
    libnames["nway"]           = "/motif_libraries_new/nway.db";
    libnames["unique_twoway"]  = "/motif_libraries_new/unique_twoway.db";
    libnames["bp_steps"]       = "/motif_libraries_new/bp_steps.db";
    
    return libnames;
}


MotifOP
MotifSqliteLibrary::get(
    String const & name,
    String const & end_id,
    String const & end_name) {
    
    String query = _generate_query(name, end_id, end_name);
    connection_.query(query);
    auto row = connection_.next();
    
    if(row->data.length() == 0) {
        throw std::runtime_error("query returned no rows");
    }
    
    if(data_.find(row->id) == data_.end() ) {
        data_[row->id] = std::make_shared<Motif>(row->data,
                                                 ResidueTypeSetManager::getInstance().residue_type_set());
    }
    
    return std::make_shared<Motif>(data_[row->id]->copy());
    
}


String
MotifSqliteLibrary::_generate_query(
    String const & name,
    String const & end_name,
    String const & end_id) {
    
    String s = "SELECT * from data_table WHERE ";
    Strings adds;
    if(name.length() > 0)     { adds.push_back("name='"+name+"' "); }
    if(end_name.length() > 0) { adds.push_back("end_name'"+end_name+"' "); }
    if(end_id.length() > 0)   { adds.push_back("end_id'"+end_id+"' "); }
    int i = 0;
    for (auto const & add : adds) {
        s += add;
        if (i != adds.size()-1) { s += "AND "; }
    }
    return s;
    
}
