//
//  motif_ensemble_sqlite_connection.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "resources/motif_ensemble_sqlite_connection.h"
#include "base/settings.h"

namespace resources {

MotifEnsembleSqliteDataOP const &MotifEnsembleSqliteConnection::next() {
  if (rc_ != SQLITE_ROW) {
    sqlite3_finalize(stmt_);
    data_->data = "";
    return data_;
  }
  data_->data =
      String(reinterpret_cast<const char *>(sqlite3_column_text(stmt_, 0)));
  data_->name =
      String(reinterpret_cast<const char *>(sqlite3_column_text(stmt_, 1)));
  data_->id =
      String(reinterpret_cast<const char *>(sqlite3_column_text(stmt_, 2)));
  rc_ = sqlite3_step(stmt_);

  return data_;
}

} // namespace resources