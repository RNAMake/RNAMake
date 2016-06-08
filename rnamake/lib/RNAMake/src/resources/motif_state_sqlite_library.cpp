//
//  motif_state_sqlite_library.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif/motif_state.h"
#include "resources/motif_state_sqlite_library.h"

StringStringMap
MotifStateSqliteLibrary::get_libnames() {
    StringStringMap libnames;
    
    libnames["ideal_helices"]  = "/motif_state_libraries/ideal_helices.db";
    libnames["ideal_helices_min"]  = "/motif_state_libraries/ideal_helices_min.db";
    libnames["twoway"]         = "/motif_state_libraries/twoway.db";
    libnames["tcontact"]       = "/motif_state_libraries/tcontact.db";
    libnames["hairpin"]        = "/motif_state_libraries/hairpin.db";
    libnames["nway"]           = "/motif_state_libraries/nway.db";
    libnames["unique_twoway"]  = "/motif_state_libraries/unique_twoway.db";
    libnames["bp_steps"]       = "/motif_state_libraries/bp_steps.db";
    
    return libnames;
}


MotifStateOP
MotifStateSqliteLibrary::get(
    String const & name,
    String const & end_id,
    String const & end_name,
    String const & id) {
    
    String query = _generate_query(name, end_id, end_name, id);
    connection_.query(query);
    auto row = connection_.next();
    
    if(row->data.length() == 0) {
        throw SqliteLibraryException(query + ": returned no rows");
    }
    
    if(data_.find(row->id) == data_.end() ) {
        data_[row->id] = std::make_shared<MotifState>(row->data);
                                                 
    }
    
    connection_.clear();
    
    return std::make_shared<MotifState>(*data_[row->id]);
    
}

MotifStateOPs
MotifStateSqliteLibrary::get_multi(
    String const & name,
    String const & end_id,
    String const & end_name,
    String const & id) {
    
    MotifStateOPs motif_states;
    String query = _generate_query(name, end_id, end_name, id);
    connection_.query(query);
    auto row = connection_.next();
    
    if(row->data.length() == 0) {
        throw SqliteLibraryException(query + ": returned no rows");
    }
    
    while(row->data.length() != 0) {
        if(data_.find(row->id) == data_.end()) {
            data_[row->id] = std::make_shared<MotifState>(row->data);

        }
        
        motif_states.push_back(std::make_shared<MotifState>(*data_[row->id]));
        row = connection_.next();
    }
    
    connection_.clear();
    
    return motif_states;
    
}

int
MotifStateSqliteLibrary::contains(
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



MotifStateOP
MotifStateSqliteLibrary::get_random() {
    int pos = rng_.randrange(max_size_);
    return get("", "", "", std::to_string(pos));
    
}

void
MotifStateSqliteLibrary::load_all(
    int limit) {
    
    int count = 0;
    for(int i = 0; i < max_size_; i++) {
        get("", "", "", std::to_string(i));
        if (count > limit) { break; }
        count++;
    }
}


String
MotifStateSqliteLibrary::_generate_query(
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