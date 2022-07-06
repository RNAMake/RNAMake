//
// Created by Joseph Yesselman on 1/5/18.
//

#include <resource_management/segment_sqlite_library.h>

namespace resource_management {

//void SegmentSqliteLibrary::get_segment(StringStringMap const &args) const {

  //_generate_query(_retrieved_columns, args);
  /*auto row = conn_.get_first_row(query_string_);
  if(row == nullptr) {
      throw SqliteLibraryException("cannot segment: " +
  base::string_map_to_string(args));
  }
  int id = row->at(0);
  auto seg = all_atom::SegmentOP(nullptr);
  if(segments_.find(id) != segments_.end() ) {
      seg = std::make_shared<all_atom::Segment>(*segments_[id]);
      seg->new_uuids();
  }
  else {
      std::vector<uint8_t> blob = row->at(1);
      auto compressed_str = String(blob.begin(), blob.end());
      auto depressed_str  = base::gzip::decompress(compressed_str.c_str(),
  compressed_str.size()); auto j = json::JSON::Load(depressed_str); seg =
  std::make_shared<all_atom::Segment>(j, rts_); segments_[id] = seg;
  }

  return seg;  */
//}

bool SegmentSqliteLibrary::contains_segment(StringStringMap const &args) const {
  /*_generate_query(_retrieved_columns, args);
  auto row = conn_.get_first_row(query_string_);
  if(row != nullptr) { return true;  }
  else               { return false; }      */
  return true;
}

} // namespace resource_management