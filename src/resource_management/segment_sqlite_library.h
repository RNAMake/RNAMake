//
// Created by Joseph Yesselman on 1/5/18.
//

#ifndef RNAMAKE_NEW_MOTIF_SQLITE_LIBRARY_H
#define RNAMAKE_NEW_MOTIF_SQLITE_LIBRARY_H

#include <resource_management/sqlite_library.h>
#include <structure/all_atom/segment.hpp>

namespace resource_management {

class SegmentSqliteLibrary : public SqliteLibrary {
public:
  SegmentSqliteLibrary(String const &db_path, String const &table_name)
      : SqliteLibrary(db_path, table_name),
        _retrieved_columns(Strings{"id", "data"}) {}

  ~SegmentSqliteLibrary() override = default;

public:
  // TODO throw if no segment exists 
  structure::all_atom::Segment get_segment(StringStringMap const & args) const {
    _generate_query(_retrieved_columns, args);
    auto row = conn_.get_first_row(query_string_);
    return structure::all_atom::get_segment_from_str(row[1].get_str());
  }

  bool contains_segment(StringStringMap const &) const;

protected:
  Strings _retrieved_columns;
  mutable std::map<int, structure::all_atom::Segment> segments_; // ack this
                                                                 // seems bad
};

} // namespace resource_management
#endif // RNAMAKE_NEW_MOTIF_SQLITE_LIBRARY_H
