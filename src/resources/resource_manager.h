//
// Created by Joseph Yesselman on 1/14/18.
//

#ifndef RNAMAKE_NEW_RESOURCE_MANAGER_H
#define RNAMAKE_NEW_RESOURCE_MANAGER_H

#include <dirent.h>

#include <structure/segment_factory.h>
#include <resources/segment_sqlite_library.h>

namespace resources {

class ResourceManagerException : public std::runtime_error {
public:
    /**
     * Standard constructor for ResourceManagerException
     * @param   message   Error message for ResourceManager
     */
    ResourceManagerException(String const & message) :
            std::runtime_error(message) {}
};

class ResourceManager {
public:
    ResourceManager():
            rts_(structure::ResidueTypeSet()),
            seg_f_(structure::SegmentFactory(rts_)) {
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
    structure::SegmentOP
    segment_from_pdb(
            String const & pdb_path,
            util::SegmentType segment_type = util::SegmentType::SEGMENT,
            bool rebuild_x3dna_files = true) const {
        return seg_f_.segment_from_pdb(pdb_path, segment_type, rebuild_x3dna_files);
    }

    inline
    structure::SegmentOPs
    all_segments_from_pdb(
            String const & pdb_path,
            util::SegmentType segment_type = util::SegmentType::SEGMENT,
            bool rebuild_x3dna_files = true) const {
        return seg_f_.all_segments_from_pdb(pdb_path, segment_type, rebuild_x3dna_files);
    }

    inline
    structure::SegmentOP
    segment_from_components(
            String const & name,
            structure::Structure const & rna_struc,
            structure::Basepairs const & basepairs,
            structure::Structure const & proteins,
            structure::Structure const & small_molecules,
            util::SegmentType segment_type = util::SegmentType::SEGMENT) const {
        return seg_f_.segment_from_components(name, rna_struc, basepairs, proteins, small_molecules, segment_type);
    }

public: // get segments

    inline
    structure::SegmentOP
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
    structure::ResidueTypeSet const &
    get_residue_type_set() { return rts_; }

private:
    structure::ResidueTypeSet rts_;
    structure::SegmentFactory seg_f_;
    std::vector<std::shared_ptr<SegmentSqliteLibrary>> sqlite_libraries_;
};

}


#endif //RNAMAKE_NEW_RESOURCE_MANAGER_H
