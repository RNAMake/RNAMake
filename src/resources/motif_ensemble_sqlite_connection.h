//
//  motif_ensemble_sqlite_connection.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_ensemble_sqlite_connection__
#define __RNAMake__motif_ensemble_sqlite_connection__

#include <stdio.h>

//RNAMake Libraries
#include "base/types.hpp"
#include "util/sqlite3_connection.h"

namespace resources {

struct MotifEnsembleSqliteData {
    MotifEnsembleSqliteData() :
            data(""), name(""), id("0") {}

    String data, name, id;

};

typedef std::shared_ptr<MotifEnsembleSqliteData> MotifEnsembleSqliteDataOP;

class MotifEnsembleSqliteConnection : public util::Sqlite3Connection {
public:
    MotifEnsembleSqliteConnection() {}

    MotifEnsembleSqliteConnection(String const & path) :
            util::Sqlite3Connection(path),
            data_(std::make_shared<MotifEnsembleSqliteData>()) {}


public:

    MotifEnsembleSqliteDataOP const &
    next();

private:
    MotifEnsembleSqliteDataOP data_;


};

}

#endif /* defined(__RNAMake__motif_ensemble_sqlite_connection__) */
