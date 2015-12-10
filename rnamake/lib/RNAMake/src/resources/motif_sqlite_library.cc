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
    String const & end_name,
    String const & id) {
    
    String query = _generate_query(name, end_id, end_name, id);
    connection_.query(query);
    auto row = connection_.next();
    
    if(row->data.length() == 0) {
        throw std::runtime_error("query returned no rows");
    }
    
    if(data_.find(row->id) == data_.end() ) {
        data_[row->id] = std::make_shared<Motif>(row->data,
                                                 ResidueTypeSetManager::getInstance().residue_type_set());
    }
    
    return std::make_shared<Motif>(*data_[row->id]);
    
}

MotifOPs
MotifSqliteLibrary::get_multi(
    String const & name,
    String const & end_id,
    String const & end_name,
    String const & id) {
    
    MotifOPs motifs;
    String query = _generate_query(name, end_id, end_name, id);
    connection_.query(query);
    auto row = connection_.next();
    
    if(row->data.length() == 0) {
        throw std::runtime_error("query returned no rows");
    }
    
    while(row->data.length() != 0) {
        if(data_.find(row->id) == data_.end()) {
            data_[row->id] = std::make_shared<Motif>(row->data,
                                                     ResidueTypeSetManager::getInstance().residue_type_set());
        }
        
        motifs.push_back(std::make_shared<Motif>(*data_[row->id]));
        row = connection_.next();
    }
    
    return motifs;
    
}

int
MotifSqliteLibrary::contains(
    String const & name,
    String const & end_id,
    String const & end_name,
    String const & id) {
    
    String query = _generate_query(name, end_id, end_name, id);
    connection_.query(query);
    auto row = connection_.next();
    
    if(row->data.length() == 0) { return 0; }
    else                        { return 1; }

    
}



MotifOP
MotifSqliteLibrary::get_random() {
    int pos = rng_.randrange(max_size_);
    return get("", "", "", std::to_string(pos));
    
}

void
MotifSqliteLibrary::load_all(
    int limit) {
    
    int count = 0;
    for(int i = 0; i < max_size_; i++) {
        get("", "", "", std::to_string(i));
        if (count > limit) { break; }
        count++;
    }
}


String
MotifSqliteLibrary::_generate_query(
    String const & name,
    String const & end_id,
    String const & end_name,
    String const & id) {
    
    String s = "SELECT * from data_table WHERE ";
    Strings adds;
    if(name.length() > 0)     { adds.push_back("name='"+name+"' "); }
    if(end_name.length() > 0) { adds.push_back("end_name='"+end_name+"' "); }
    if(end_id.length() > 0)   { adds.push_back("end_id='"+end_id+"' "); }
    if(id.length() > 0)       { adds.push_back("id='"+id+"' "); }

    int i = 0;
    for (auto const & add : adds) {
        s += add;
        if (i != adds.size()-1) { s += "AND "; }
        i++;
    }
    return s;
    
}
