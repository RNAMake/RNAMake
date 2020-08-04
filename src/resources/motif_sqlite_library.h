//
//  motif_sqlite_library.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_sqlite_library__
#define __RNAMake__motif_sqlite_library__

#include <stdio.h>
#include <iostream>

#include "util/random_number_generator.h"
#include "motif/motif.h"
#include "resources/motif_sqlite_connection.h"
#include "resources/sqlite_library.h"
#include "base/settings.h"

namespace resources {

static String dummy_name = "", dummy_end_id = "", dummy_end_name = "", dummy_id = "";

class MotifSqliteLibrary : public SqliteLibrary {
public:

    MotifSqliteLibrary(
            String const & libname) {

        libnames_ = get_libnames();
        rng_ = util::RandomNumberGenerator();
        auto path = _get_path(libname);
        MotifSqliteConnection conn(path);
        connection_ = conn;
        max_size_ = connection_.count();
        //max_size_ = 1;
    }

    ~MotifSqliteLibrary() {}

public: //iterator stuff

    class iterator {
    public:
        iterator(
                std::map<String, motif::MotifOP>::iterator const & i) :
                i_(i) {}

        iterator operator++() {
            i_++;
            return *this;
        }

        motif::MotifOP const & operator*() { return i_->second; }

        bool operator==(iterator const & rhs) const { return i_ == rhs.i_; }

        bool operator!=(iterator const & rhs) const { return i_ != rhs.i_; }

    private:
        std::map<String, motif::MotifOP>::iterator i_;

    };


    iterator begin() { return iterator(data_.begin()); }

    iterator end() { return iterator(data_.end()); }


public:

    static
    StringStringMap
    get_libnames();

    motif::MotifOP
    get(
            String const & name = dummy_name,
            String const & end_id = dummy_end_id,
            String const & end_name = dummy_name,
            String const & id = dummy_id);

    motif::MotifOPs
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

    motif::MotifOP
    get_random();

    void
    load_all(
            int limit = 99999);


private:

    String
    _generate_query(
            String const &,
            String const &,
            String const &,
            String const &);


private:

    MotifSqliteConnection connection_;
    std::map<String, motif::MotifOP> data_;
    util::RandomNumberGenerator rng_;

};


typedef std::shared_ptr<MotifSqliteLibrary> MotifSqliteLibraryOP;

}

#endif /* defined(__RNAMake__motif_sqlite_library__) */
