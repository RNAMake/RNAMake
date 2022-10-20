//
//  sqlite3_connection.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <base/paths.hpp>
#include <resources/segment_sqlite_connection.h>

namespace resource_management {

SegmentSqliteDataOP const &SegmentSqliteConnection::next() {
  if (_rc != SQLITE_ROW) {
    sqlite3_finalize(_stmt);
    _data->data = "";
    return _data;
  }
  _data->data =
      String(reinterpret_cast<const char *>(sqlite3_column_text(_stmt, 0)));
  _data->name =
      String(reinterpret_cast<const char *>(sqlite3_column_text(_stmt, 1)));
  _data->end_name =
      String(reinterpret_cast<const char *>(sqlite3_column_text(_stmt, 2)));
  _data->end_id =
      String(reinterpret_cast<const char *>(sqlite3_column_text(_stmt, 3)));
  _data->id =
      String(reinterpret_cast<const char *>(sqlite3_column_text(_stmt, 4)));
  _rc = sqlite3_step(_stmt);
  return _data;
}

SegmentSqliteDataOP const &SegmentSqliteConnection::contains() {
  if (_rc != SQLITE_ROW) {
    sqlite3_finalize(_stmt);
    _data->data = "";
    return _data;
  }
  _data->data =
      String(reinterpret_cast<const char *>(sqlite3_column_text(_stmt, 0)));
  _data->name =
      String(reinterpret_cast<const char *>(sqlite3_column_text(_stmt, 1)));
  _data->end_name =
      String(reinterpret_cast<const char *>(sqlite3_column_text(_stmt, 2)));
  _data->end_id =
      String(reinterpret_cast<const char *>(sqlite3_column_text(_stmt, 3)));
  _data->id =
      String(reinterpret_cast<const char *>(sqlite3_column_text(_stmt, 4)));

  sqlite3_finalize(_stmt);
  return _data;
}

} // namespace resource_management
