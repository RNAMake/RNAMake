//
// Created by Joseph Yesselman on 12/30/17.
//

#ifndef RNAMAKE_NEW_SQLITE_DATABASE_H
#define RNAMAKE_NEW_SQLITE_DATABASE_H

#include <sqlite3.h>
#include <base/types.h>
#include <base/file_io.h>

namespace util {
namespace sqlite {

class Database {
public:
    Database(
            String name) :
            name_(name),
            db_(nullptr),
            created_(false){ open(); }

    ~Database() {
        // close the db
        (void) close();
    }

public:
    // close the database
    int close() {
        int err = sqlite3_close(db_);
        open_ = false;
        db_ = nullptr;
        return err;
    }

    // returns true if the database is open
    inline
    bool
    is_open() const { return open_; }

    inline
    bool
    is_created() const { return created_; }

    // SQLite3 access
    inline
    sqlite3
    *get() const { return db_; }

    inline
    sqlite3
    *operator()() const { return db_; }

    inline
    String const &
    get_name() const { return name_; }


private:
    // open (connect) the database
    int open(int flags = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE) {
        if(!base::file_exists(name_)) { created_ = true; }

        int err = sqlite3_open_v2(name_.c_str(), &db_, flags, nullptr);
        open_ = !err;
        return err;
    }

private:
    Database & operator=(const Database &) = delete;

private:
    sqlite3 *db_;    // associated db
    String const name_;  // db filename
    bool open_;  // db open status
    bool created_;
};

}
}

#endif //RNAMAKE_NEW_SQLITE_DATABASE_H
