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
#include "resources/motif_sqlite_library.h"

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
    
    void
    add_motif(
        String const &);
    
    
protected:
    ResourceManager() { //Prevent construction
        mf_ = MotifFactory();
        mlibs_ = std::map<String, MotifSqliteLibraryOP>();
        
        for(auto const & kv : MotifSqliteLibrary::get_libnames()) {
            mlibs_[kv.first] = std::make_shared<MotifSqliteLibrary>(kv.first);
        }
    
    }
    
    ResourceManager(ResourceManager const &); //Prevent construction
    void operator= (ResourceManager const &);
    
private:
    ~ResourceManager() {}
    
private:
    std::map<String, MotifSqliteLibraryOP> mlibs_;
    MotifFactory mf_;

    
};


#endif /* defined(__RNAMake__resource_manager__) */
