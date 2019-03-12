//
//  motif_state_ensemble_sqlite_library.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_ensemble_sqlite_library__
#define __RNAMake__motif_state_ensemble_sqlite_library__

#include <stdio.h>


#include <stdio.h>

#include "util/random_number_generator.h"
#include "motif/motif_state_ensemble.h"
#include "resources/motif_ensemble_sqlite_connection.h"
#include "resources/motif_sqlite_library.h"

class MotifStateEnsembleSqliteLibrary : public SqliteLibrary {
public:
    
    MotifStateEnsembleSqliteLibrary(
        String const & libname) {
        
        libnames_ = get_libnames();
        rng_ = util::RandomNumberGenerator();
        auto path = _get_path(libname);
        MotifEnsembleSqliteConnection conn(path);
        connection_ = conn;
        max_size_ = connection_.count();        
    }
    
    ~MotifStateEnsembleSqliteLibrary() {}
    
public: //iterator stuff
    
    class iterator {
    public:
        iterator(
            std::map<String, MotifStateEnsembleOP>::iterator const & i ):
        i_(i)
        {}
        
        iterator operator++() { ++i_; return *this; }
        MotifStateEnsembleOP const & operator*() { return i_->second; }
        bool operator== (iterator const & rhs) const { return i_ == rhs.i_; }
        bool operator!= (iterator const & rhs) const { return i_ != rhs.i_; }
        
    private:
        std::map<String, MotifStateEnsembleOP>::iterator i_;
        
    };
    
    
    iterator begin() { return iterator(data_.begin()); }
    iterator end()   { return iterator(data_.end()); }
    
    
public:
    
    static
    StringStringMap
    get_libnames();
    
    MotifStateEnsembleOP
    get(
        String const & name = dummy_name,
        String const & id = dummy_id);
    
    MotifStateEnsembleOPs
    get_multi(
        String const & name = dummy_name,
        String const & id = dummy_id);
    
    int
    contains(
        String const & name = dummy_name,
        String const & id = dummy_id);
    
    MotifStateEnsembleOP
    get_random();
    
    void
    load_all(
        int limit=99999);
    
    
private:
    
    String
    _generate_query(
        String const &,
        String const &);
    
private:
    
    MotifEnsembleSqliteConnection connection_;
    std::map<String, MotifStateEnsembleOP> data_;
    util::RandomNumberGenerator rng_;
    
};


typedef std::shared_ptr<MotifStateEnsembleSqliteLibrary> MotifStateEnsembleSqliteLibraryOP;


#endif /* defined(__RNAMake__motif_state_ensemble_sqlite_library__) */
