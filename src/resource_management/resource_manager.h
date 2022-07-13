//
// Created by Joseph Yesselman on 1/14/18.
//

#ifndef RNAMAKE_NEW_RESOURCE_MANAGER_H
#define RNAMAKE_NEW_RESOURCE_MANAGER_H

#include <filesystem>

#include <resource_management/exception.hpp>
#include <resource_management/segment_sqlite_library.h>

namespace resource_management {

class ResourceManager {
public:
  ResourceManager() {
    auto path = base::path::resources_path() + "/motif_libraries_new/";
    if (!std::filesystem::is_directory(path)) {
      String msg = "cannot locate sqlite3 resource directory! Something is "
                   "very wrong!";
      throw base::ResourceException(msg);
    }
    for (const auto &entry : std::filesystem::directory_iterator(path)) {
      auto fname = String(entry.path());
      LOG_DEBUG << "resource manager loaded: " << fname;
      _sqlite_libraries.emplace_back(fname, "data_table");
    }
  }

  ~ResourceManager() = default;
  // only want one resource manager active at a time
  ResourceManager(const ResourceManager &) = delete;
  // only want one resrouce manager active at a time
  ResourceManager &operator=(const ResourceManager &) { return *this; }

public: // load new segments from pdbs and components
  /*inline all_atom::SegmentOP
  segment_from_pdb(String const &pdb_path,
                   util::SegmentType segment_type = util::SegmentType::SEGMENT,
                   bool rebuild_x3dna_files = true) const {
    return seg_f_.segment_from_pdb(pdb_path, segment_type, rebuild_x3dna_files);
  }

  inline all_atom::SegmentOPs all_segments_from_pdb(
      String const &pdb_path,
      util::SegmentType segment_type = util::SegmentType::SEGMENT,
      bool rebuild_x3dna_files = true) const {
    return seg_f_.all_segments_from_pdb(pdb_path, segment_type,
                                        rebuild_x3dna_files);
  }

  inline all_atom::SegmentOP segment_from_components(
      String const &name, all_atom::Structure const &rna_struc,
      all_atom::Basepairs const &basepairs, all_atom::Structure const &proteins,
      all_atom::Structure const &small_molecules,
      util::SegmentType segment_type = util::SegmentType::SEGMENT) const {
    return seg_f_.segment_from_components(name, rna_struc, basepairs, proteins,
                                          small_molecules, segment_type);
  } */

public: // get segments
  structure::all_atom::SegmentOP get_segment(const SegmentInfo &seg_info) {
    for (const auto &seg_lib : _sqlite_libraries) {
      if (seg_lib.contains_segment(seg_info)) {
        return std::make_shared<structure::all_atom::Segment>(
            seg_lib.get_segment(seg_info));
      }
    }
    throw ResourceManagementException("cannot find segment: " +
                                      seg_info.get_str());
  }

  [[nodiscard]] inline bool
  contains_segment(const SegmentInfo &seg_info) const {
    for (auto &seg_lib : _sqlite_libraries) {
      if (seg_lib.contains_segment(seg_info)) {
        return true;
      }
    }
    return false;
  }

private:
  // all_atom::ResidueTypeSet rts_;
  // all_atom::SegmentFactory seg_f_;
  std::vector<SegmentSqliteLibrary> _sqlite_libraries;
};

} // namespace resource_management

#endif // RNAMAKE_NEW_RESOURCE_MANAGER_H
