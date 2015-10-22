//
//  motif_state_ensemble_sqlite_library.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_state_ensemble_sqlite_library.h"

StringStringMap
MotifStateEnsembleSqliteLibrary::get_libnames() {
    StringStringMap libnames;
    
    libnames["all_bp_steps"]   = "/motif_state_ensemble_libraries/all_bp_steps.db";
    libnames["bp_steps"]       = "/motif_state_ensemble_libraries/bp_steps.db";
    libnames["twoway"]         = "/motif_state_ensemble_libraries/twoway.db";
    libnames["tcontact"]       = "/motif_state_ensemble_libraries/tcontact.db";
    libnames["hairpin"]        = "/motif_state_ensemble_libraries/hairpin.db";
    libnames["nway"]           = "/motif_state_ensemble_libraries/nway.db";
    
    return libnames;
}


MotifStateEnsembleOP
MotifStateEnsembleSqliteLibrary::get(
    String const & name,
    String const & id) {
    
    String query = _generate_query(name, id);
    connection_.query(query);
    auto row = connection_.next();
    
    if(row->data.length() == 0) {
        throw std::runtime_error("query returned no rows");
    }
    
    if(data_.find(row->id) == data_.end() ) {
        data_[row->id] = std::make_shared<MotifStateEnsemble>(row->data);
        
    }
    
    return std::make_shared<MotifStateEnsemble>(data_[row->id]->copy());
    
}

MotifStateEnsembleOPs
MotifStateEnsembleSqliteLibrary::get_multi(
    String const & name,
    String const & id) {
    
    MotifStateEnsembleOPs motif_state_ensembles;
    String query = _generate_query(name, id);
    connection_.query(query);
    auto row = connection_.next();
    
    if(row->data.length() == 0) {
        throw std::runtime_error("query returned no rows");
    }
    
    while(row->data.length() != 0) {
        if(data_.find(row->id) == data_.end()) {
            data_[row->id] = std::make_shared<MotifStateEnsemble>(row->data);
            
        }
        
        motif_state_ensembles.push_back(std::make_shared<MotifStateEnsemble>(data_[row->id]->copy()));
        row = connection_.next();
    }
    
    return motif_state_ensembles;
    
}

int
MotifStateEnsembleSqliteLibrary::contains(
    String const & name,
    String const & id) {
    
    String query = _generate_query(name, id);
    connection_.query(query);
    auto row = connection_.next();
    
    if(row->data.length() == 0) { return 0; }
    else                        { return 1; }
    
    
}

MotifStateEnsembleOP
MotifStateEnsembleSqliteLibrary::get_random() {
    int pos = rng_.randrange(max_size_);
    return get("", std::to_string(pos));
    
}

void
MotifStateEnsembleSqliteLibrary::load_all(
    int limit) {
    
    int count = 0;
    for(int i = 0; i < max_size_; i++) {
        get("", std::to_string(i));
        if (count > limit) { break; }
        count++;
    }
}

String
MotifStateEnsembleSqliteLibrary::_generate_query(
    String const & name,
    String const & id) {
    
    String s = "SELECT * from data_table WHERE ";
    Strings adds;
    if(name.length() > 0)     { adds.push_back("name='"+name+"' "); }
    if(id.length() > 0)       { adds.push_back("id='"+id+"' "); }
    
    int i = 0;
    for (auto const & add : adds) {
        s += add;
        if (i != adds.size()-1) { s += "AND "; }
        i++;
    }
    return s;
    
}