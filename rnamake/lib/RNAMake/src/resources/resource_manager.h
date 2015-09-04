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


class ResourceManager {
public:
    static ResourceManager & getInstance() {
        static ResourceManager instance;
        return instance;
    }
    
public:
    
    MotifOP
    get_motif(
        String const & name = dummy_name,
        String const & end_id = dummy_end_id,
        String const & end_name = dummy_name,
        String const & id = dummy_id);
    
    MotifStateOP
    get_state(
        String const & name = dummy_name,
        String const & end_id = dummy_end_id,
        String const & end_name = dummy_name,
        String const & id = dummy_id);
    
    MotifStateEnsembleOP
    get_motif_state_ensemble(
        String const & name = dummy_name,
        String const & id = dummy_id);
    
    void
    add_motif(
        String const &);
    
    
protected:
    ResourceManager() { //Prevent construction
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
    
    ResourceManager(ResourceManager const &); //Prevent construction
    void operator= (ResourceManager const &);
    
private:
    ~ResourceManager() {}
    
private:
    std::map<String, MotifSqliteLibraryOP> mlibs_;
    std::map<String, MotifStateSqliteLibraryOP> ms_libs_;
    std::map<String, MotifStateEnsembleSqliteLibraryOP> mse_libs_;
    MotifFactory mf_;
    AddedMotifLibrary added_motifs_;

    
};


#endif /* defined(__RNAMake__resource_manager__) */
