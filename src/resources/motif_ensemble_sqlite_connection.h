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

// RNAMake Libraries
#include "base/types.hpp"
#include "util/sqlite/connection.hpp"

namespace resources {

struct MotifEnsembleSqliteData {
  MotifEnsembleSqliteData() : data(""), name(""), id("0") {}

  String data, name, id;
};

typedef std::shared_ptr<MotifEnsembleSqliteData> MotifEnsembleSqliteDataOP;

class MotifEnsembleSqliteConnection : public util::sqlite::Connection {
public:
  MotifEnsembleSqliteConnection() {}

  MotifEnsembleSqliteConnection(String const &path)
      : util::sqlite::Connection(path),
        _data(std::make_shared<MotifEnsembleSqliteData>()) {}

public:
  MotifEnsembleSqliteDataOP const &next();

private:
  MotifEnsembleSqliteDataOP _data;
};

} // namespace resources

#endif /* defined(__RNAMake__motif_ensemble_sqlite_connection__) */
