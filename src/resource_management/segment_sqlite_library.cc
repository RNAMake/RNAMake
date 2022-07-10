//
// Created by Joseph Yesselman on 1/5/18.
//

#include <resource_management/segment_sqlite_library.h>

namespace resource_management {

// TODO throw if no segment exists
structure::all_atom::Segment
SegmentSqliteLibrary::get_segment(const SegmentInfo & seg_info) const {
  _generate_query(_retrieved_columns, seg_info.get_dict());
  try {
    auto row = _conn.get_first_row(_query_string);
    return structure::all_atom::get_segment_from_str(row[1].get_str());
  }
  catch(const util::SqliteException & e) {
    throw ResourceManagementException(e.what());
  }
}

bool SegmentSqliteLibrary::contains_segment(const SegmentInfo & seg_info) const {
  if(_does_query_return_rows(seg_info.get_dict())) {
    return true;
  }
  else {
    return false;
  }
}

} // namespace resource_management