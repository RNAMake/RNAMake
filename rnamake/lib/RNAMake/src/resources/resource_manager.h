//
//  resource_manager.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__resource_manager__
#define __RNAMake__resource_manager__

#include <stdio.h>

//RNAMake Headers
#include "base/types.h"
#include "motif/motif_factory.h"
#include "motif/motif_state.h"
#include "resources/motif_sqlite_library.h"
#include "resources/motif_state_sqlite_library.h"
#include "resources/motif_state_ensemble_sqlite_library.h"
#include "resources/added_motif_library.h"


class ResourceManagerException : public std::runtime_error {
public:
    ResourceManagerException(
        String const & message) :
    std::runtime_error(message)
    {}
};


class RM { //RM for ResourceManager
public:
    static RM & instance() {
        static RM instance;
        return instance;
    }
    
public:
    
    MotifOP
    motif(
        String const & name = dummy_name,
        String const & end_id = dummy_end_id,
        String const & end_name = dummy_name);
    
    MotifStateOP
    motif_state(
        String const & name = dummy_name,
        String const & end_id = dummy_end_id,
        String const & end_name = dummy_name);
    
    MotifStateEnsembleOP
    motif_state_ensemble(
        String const & name = dummy_name);
    
    void
    add_motif(
        String const &);
    
    
protected:
    RM() { //Prevent construction
        mf_ = MotifFactory();
        added_motifs_ = AddedMotifLibrary();
        mlibs_ = std::map<String, MotifSqliteLibraryOP>();
        ms_libs_ = std::map<String, MotifStateSqliteLibraryOP>();
        mse_libs_ = std::map<String, MotifStateEnsembleSqliteLibraryOP>();
        
        for(auto const & kv : MotifSqliteLibrary::get_libnames()) {
            mlibs_[kv.first] = std::make_shared<MotifSqliteLibrary>(kv.first);
        }
        
        for(auto const & kv : MotifStateSqliteLibrary::get_libnames()) {
            ms_libs_[kv.first] = std::make_shared<MotifStateSqliteLibrary>(kv.first);
        }
        
        for(auto const & kv : MotifStateEnsembleSqliteLibrary::get_libnames()) {
            mse_libs_[kv.first] = std::make_shared<MotifStateEnsembleSqliteLibrary>(kv.first);
        }
    
    }
    
    RM(RM const &); //Prevent construction
    void operator= (RM const &);
    
private:
    ~RM() {}
    
private:
    std::map<String, MotifSqliteLibraryOP> mlibs_;
    std::map<String, MotifStateSqliteLibraryOP> ms_libs_;
    std::map<String, MotifStateEnsembleSqliteLibraryOP> mse_libs_;
    MotifFactory mf_;
    AddedMotifLibrary added_motifs_;

    
};

inline
MotifOP
get_motif_from_resource_manager(
    String const & name = dummy_name,
    String const & end_id = dummy_end_id,
    String const & end_name = dummy_name) {
    
    return RM::instance().motif(name, end_id, end_name);
    
}


#endif /* defined(__RNAMake__resource_manager__) */
