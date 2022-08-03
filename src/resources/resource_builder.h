//
// Created by Joseph Yesselman on 1/4/18.
//

#ifndef RNAMAKE_NEW_RESOURCE_BUILDER_H
#define RNAMAKE_NEW_RESOURCE_BUILDER_H

#include <base/settings.h>
#include <structure/segment_factory.h>
#include <util/sqlite/connection.h>
#include <base/types.hpp>
#include <util/sqlite/connection.hpp>

namespace resources {

class ResourceBuilderException : public std::runtime_error {
public:
  /**
   * Standard constructor for ResourceBuilderException
   * @param   message   Error message for SqliteLibrary
   */
  ResourceBuilderException(String const &message)
      : std::runtime_error(message) {}
};

class ResourceBuilder {
public:
  inline ResourceBuilder(structure::SegmentFactory &seg_f,
                         String const &motif_dirs_path)
      : _seg_f(seg_f), _motif_dirs_path(motif_dirs_path),
        _motif_table(_generate_motif_table_details()) {
    _start_insert_str =
        "INSERT INTO data_table (id, data, name, end_name, end_id) VALUES (";
  }

public:
  void build_ideal_helices();

  void build_basic_libraries();

  // private:
  //     void
  //     _insert_segment_to_motif_table(
  //             structure::Segment const &,
  //             int,
  //             util::sqlite::Connection &);

  util::sqlite::TableDetails _generate_motif_table_details() {
    auto motif_table = util::sqlite::TableDetails("data_table");
    motif_table.add_column("data", "BLOB");
    motif_table.add_column("name", "TEXT");
    motif_table.add_column("end_name", "TEXT");
    motif_table.add_column("end_id", "TEXT");
    motif_table.add_column("id", "INT", true);
    return motif_table;
  }

  void _build_basic_library(String const &, util::SegmentType,
                            std::map<String, int> const &);

private:
  structure::SegmentFactory &_seg_f;
  String _motif_dirs_path;
  util::sqlite::TableDetails _motif_table;
  std::vector<uint8_t> _blob;
  String _compressed_str;
  String _start_insert_str;
};

} // namespace resources

#endif // RNAMAKE_NEW_RESOURCE_BUILDER_H
