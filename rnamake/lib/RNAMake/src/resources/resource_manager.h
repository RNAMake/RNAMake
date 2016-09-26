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
#include "motif/motif_ensemble.h"
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
protected:
    
    RM();
    
    RM(RM const &); //Prevent construction
    void operator= (RM const &);
    
private:

    ~RM() {}
    
public:
    
    static RM & instance() {
        static RM instance;
        return instance;
    }

public: // getting functions
    
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
    
public: // adding functions

    RNAStructureOP
    get_structure(
        String const & path,
        String name = "");
    
    void
    add_motif(
        String const & path,
        String name = "");
    
    void
    register_motif(
        MotifOP const &);
    
    void
    register_extra_motif_ensembles(
        String const &);
    
    int
    has_supplied_motif_ensemble(
        String const &,
        String const &);
    
    MotifEnsembleOP const & 
    get_supplied_motif_ensemble(
        String const &,
        String const &);
    
    
private:
    std::map<String, MotifSqliteLibraryOP> mlibs_;
    std::map<String, MotifStateSqliteLibraryOP> ms_libs_;
    std::map<String, MotifStateEnsembleSqliteLibraryOP> mse_libs_;
    std::map<String, MotifEnsembleOP> extra_me_;
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
