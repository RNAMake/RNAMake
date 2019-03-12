//
//  motif_state_sqlite_library.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_sqlite_library__
#define __RNAMake__motif_state_sqlite_library__

#include <stdio.h>

#include "util/random_number_generator.h"
#include "motif/motif_state.h"
#include "resources/motif_sqlite_connection.h"
#include "resources/motif_sqlite_library.h"

class MotifStateSqliteLibrary : public SqliteLibrary {
public:
    
    MotifStateSqliteLibrary(
        String const & libname) {
        
        libnames_ = get_libnames();
        rng_ = util::RandomNumberGenerator();
        auto path = _get_path(libname);
        MotifSqliteConnection conn(path);
        connection_ = conn;
        max_size_ = connection_.count();
        
    }
    
    ~MotifStateSqliteLibrary() {}
    
public: //iterator stuff
    
    class iterator {
    public:
        iterator(
            std::map<String, MotifStateOP>::iterator const & i ):
        i_(i)
        {}
        
        iterator operator++() { i_++; return *this; }
        MotifStateOP const & operator*() { return i_->second; }
        bool operator== (iterator const & rhs) const { return i_ == rhs.i_; }
        bool operator!= (iterator const & rhs) const { return i_ != rhs.i_; }
        
    private:
        std::map<String, MotifStateOP>::iterator i_;
        
    };
    
    
    iterator begin() { return iterator(data_.begin()); }
    iterator end()   { return iterator(data_.end()); }
    
    
public:
    
    static
    StringStringMap
    get_libnames();
    
    MotifStateOP
    get(
        String const & name = dummy_name,
        String const & end_id = dummy_end_id,
        String const & end_name = dummy_name,
        String const & id = dummy_id);
    
    MotifStateOPs
    get_multi(
        String const & name = dummy_name,
        String const & end_id = dummy_end_id,
        String const & end_name = dummy_name,
        String const & id = dummy_id);
    
    int
    contains(
        String const & name = dummy_name,
        String const & end_id = dummy_end_id,
        String const & end_name = dummy_name,
        String const & id = dummy_id);
    
    MotifStateOP
    get_random();
    
    void
    load_all(
        int limit=99999);
    
    
private:
    
    String
    _generate_query(
        String const &,
        String const &,
        String const &,
        String const &);    
    
private:
    
    MotifSqliteConnection connection_;
    std::map<String, MotifStateOP> data_;
    util::RandomNumberGenerator rng_;
    
};


typedef std::shared_ptr<MotifStateSqliteLibrary> MotifStateSqliteLibraryOP;


#endif /* defined(__RNAMake__motif_state_sqlite_library__) */
