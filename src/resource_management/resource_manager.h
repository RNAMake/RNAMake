//
// Created by Joseph Yesselman on 1/14/18.
//

#ifndef RNAMAKE_NEW_RESOURCE_MANAGER_H
#define RNAMAKE_NEW_RESOURCE_MANAGER_H

#include <dirent.h>

#include <resource_management/segment_sqlite_library.h>

/*
namespace resource_management {

class ResourceManager {
public:
    ResourceManager():
            rts_(all_atom::ResidueTypeSet()),
            seg_f_(all_atom::SegmentFactory(rts_)) {
        auto path = base::resources_path() + "/motif_libraries/";

        DIR *pDIR = opendir(path.c_str());
        struct dirent *entry;
        while ((entry = readdir(pDIR)) != NULL) {
            auto fname = String(entry->d_name);
            if(fname[0] == '.') { continue; } // remove . and ..
            auto seg_lib = std::make_shared<SegmentSqliteLibrary>(path+fname, "data_table", rts_);
            sqlite_libraries_.push_back(seg_lib);
        }

        closedir(pDIR);
        delete entry;

    }

    ~ResourceManager() {}

public:
    // only want one resource manager active at a time
    ResourceManager(const ResourceManager&) = delete;
    ResourceManager & operator = (const ResourceManager & ) { return *this; }

public:    // load new segments from pdbs and components
    inline
    all_atom::SegmentOP
    segment_from_pdb(
            String const & pdb_path,
            util::SegmentType segment_type = util::SegmentType::SEGMENT,
            bool rebuild_x3dna_files = true) const {
        return seg_f_.segment_from_pdb(pdb_path, segment_type, rebuild_x3dna_files);
    }

    inline
    all_atom::SegmentOPs
    all_segments_from_pdb(
            String const & pdb_path,
            util::SegmentType segment_type = util::SegmentType::SEGMENT,
            bool rebuild_x3dna_files = true) const {
        return seg_f_.all_segments_from_pdb(pdb_path, segment_type, rebuild_x3dna_files);
    }

    inline
    all_atom::SegmentOP
    segment_from_components(
            String const & name,
            all_atom::Structure const & rna_struc,
            all_atom::Basepairs const & basepairs,
            all_atom::Structure const & proteins,
            all_atom::Structure const & small_molecules,
            util::SegmentType segment_type = util::SegmentType::SEGMENT) const {
        return seg_f_.segment_from_components(name, rna_struc, basepairs, proteins, small_molecules, segment_type);
    }

public: // get segments

    inline
    all_atom::SegmentOP
    get_segment(
            StringStringMap const & args) const {
        for(auto & seg_lib : sqlite_libraries_) {
            if(seg_lib->contains_segment(args)) { return seg_lib->get_segment(args); }
        }

        throw ResourceManagerException("cannot find segment: " + base::string_map_to_string(args));
    }

    inline
    bool
    contains_segment(
            StringStringMap const & args) const {
        for (auto & seg_lib : sqlite_libraries_) {
            if(seg_lib->contains_segment(args)) { return true; }
        }
        return false;
    }


public:
    all_atom::ResidueTypeSet const &
    get_residue_type_set() { return rts_; }

private:
    all_atom::ResidueTypeSet rts_;
    all_atom::SegmentFactory seg_f_;
    std::vector<std::shared_ptr<SegmentSqliteLibrary>> sqlite_libraries_;
};

}
  */

#endif //RNAMAKE_NEW_RESOURCE_MANAGER_H
