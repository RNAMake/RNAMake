//
// Created by Joseph Yesselman on 1/5/18.
//

#ifndef RNAMAKE_NEW_MOTIF_SQLITE_LIBRARY_H
#define RNAMAKE_NEW_MOTIF_SQLITE_LIBRARY_H

#include <resource_management/sqlite_library.h>
#include <structure/all_atom/segment.hpp>

namespace resource_management {

struct SegmentInfo {
  String name;
  String end_name;
  String end_id;

  [[nodiscard]] bool is_valid() const {
    if (name.empty() && end_name.empty() && end_id.empty()) {
      return false;
    } else {
      return true;
    }
  }

  [[nodiscard]] StringStringMap get_dict() const {
    StringStringMap d;
    if (!name.empty()) {
      d["name"] = name;
    }
    if (!end_name.empty()) {
      d["end_name"] = end_name;
    }
    if (!end_id.empty()) {
      d["end_id"] = end_id;
    }
    return d;
  }

  [[nodiscard]] String get_str() const  {
    String s;
    if (!name.empty()) {
      s += "name: " + name;
    }
    if (!end_name.empty()) {
      if(!s.empty()) { s += " "; }
      s += "end_name: " + end_name;
    }
    if (!end_id.empty()) {
      if(!s.empty()) { s += " "; }
      s += "end_id: " + end_id;
    }
    return s;
  }
};

class SegmentSqliteLibrary : public SqliteLibrary {
public:
  SegmentSqliteLibrary(const String &db_path, const String &table_name)
      : SqliteLibrary(db_path, table_name),
        _retrieved_columns(Strings{"id", "data"}) {
  }

  //SegmentSqliteLibrary(SegmentSqliteLibrary && other) {

  //}

  ~SegmentSqliteLibrary() override = default;

public:
  structure::all_atom::Segment get_segment(const SegmentInfo &) const;

  bool contains_segment(const SegmentInfo &) const;

protected:
  Strings _retrieved_columns;
  mutable std::map<int, structure::all_atom::Segment> _segments; // ack this
                                                                 // seems bad
};

} // namespace resource_management
#endif // RNAMAKE_NEW_MOTIF_SQLITE_LIBRARY_H
