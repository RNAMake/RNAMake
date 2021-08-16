//
// Created by Joseph Yesselman on 1/4/18.
//
#include <dirent.h>

#include <base/string.h>
#include <base/gzip/compress.hpp>
#include <resources/resource_builder.h>

namespace resources {

//void
//ResourceBuilder::build_ideal_helices() {
//    auto helix_path = motif_dirs_path_ + "/helices/";
//
//    DIR *pDIR;
//    struct dirent *entry;
//
//    pDIR = opendir(helix_path.c_str());
//    // directory does not exist
//    if(ENOENT == errno) { throw ResourceBuilderException("motif helix path does not exist!"); }
//
//    auto lib_path = base::resources_path() + "/motif_libraries/ideal_helices.db";
//
//    LOGI << "building ideal helix segement library at: " + lib_path;
//
//    try { std::remove(lib_path.c_str()); }
//    catch(...) { }
//
//    auto db = util::sqlite::Database(lib_path);
//    auto conn = util::sqlite::Connection(db);
//    util::sqlite::create_table(conn, motif_table_);
//
//    auto i = 0;
//    auto insert_str = String("");
//
//    while ((entry = readdir(pDIR)) != NULL) {
//        auto fname = String(entry->d_name);
//        if(fname.length() < 10) { continue;}
//        auto spl = base::split_str_by_delimiter(fname, ".");
//        if(spl[1] != "IDEAL") { continue; }
//
//        auto path = helix_path + fname + "/" + fname + ".pdb";
//        auto segs = seg_f_.all_segments_from_pdb(path, util::SegmentType::HELIX, false);
//
//        seg_f_.align_segment_to_ref_frame(*segs[0]);
//        _insert_segment_to_motif_table(*segs[0], i, conn);
//        i++;
//    }
//
//    closedir(pDIR);
//    delete entry;
//}

void
ResourceBuilder::build_basic_libraries() {
    auto paths = std::map<String, util::SegmentType>{
            { motif_dirs_path_ + "/two_ways/" , util::SegmentType::TWOWAY_JUNCTION },
            { motif_dirs_path_ + "/hairpins/" , util::SegmentType::HAIRPIN },
            { motif_dirs_path_ + "/helices/"  , util::SegmentType::HELIX },
            { motif_dirs_path_ + "/junctions/", util::SegmentType::NWAY_JUNCTION },
            { motif_dirs_path_ + "/tertiary_contact_hairpin_hairpin/",   util::SegmentType::TC_HAIRPIN_HAIRPIN },
            { motif_dirs_path_ + "/tertiary_contact_junction_hairpin/",  util::SegmentType::TC_JUNCTION_HAIRPIN },
            { motif_dirs_path_ + "/tertiary_contact_junction_junction/", util::SegmentType::TC_JUNCTION_JUNCTION},
            { motif_dirs_path_ + "/basepair_steps/", util::SegmentType::HELIX }
    };

    auto excluded_motifs = std::map<String, int>();
    // read in files exlcuding motifs
    // motifs that I have investigated why they cannot be loaded into RNAMake
    auto path = base::resources_path() + "/motif_info/explained_failures.csv";
    auto lines = base::get_lines_from_file(path);
    for(auto const & l : lines) {
        if(l.size() < 10) { continue; }
        auto spl = base::split_str_by_delimiter(l, ",");
        excluded_motifs[spl[0]] = 1;
    }

    // motifs that have extra motifs
    path = base::resources_path() + "/motif_info/motifs_with_tcs.csv";
    lines = base::get_lines_from_file(path);
    for(auto const & l : lines) {
        if(l.size() < 10) { continue; }
        auto spl = base::split_str_by_delimiter(l, ",");
        if(spl.size() < 4) { continue; }
        excluded_motifs[spl[0]] = 1;
    }

    for(auto const & kv : paths) {
        _build_basic_library(kv.first, kv.second, excluded_motifs);
    }

}


//void
//ResourceBuilder::_insert_segment_to_motif_table(
//        structure::Segment const & seg,
//        int id,
//        util::sqlite::Connection & conn) {
//
//    auto insert_str = String("");
//    //TODO Need to change this to not use json
//    auto j = seg.get_json();
//    auto j_str = j.dump_minimized();
//    compressed_str_ = base::gzip::compress(j_str.c_str(), j_str.size());
//    blob_ = std::vector<uint8_t>(compressed_str_.begin(), compressed_str_.end());
//
//    conn.bind(1, blob_);
//    insert_str  = start_insert_str_ + base::quoted_string(std::to_string(id)) + ", ?, ";
//    insert_str += base::quoted_string(seg.get_name_str()) + ", ";
//    insert_str += base::quoted_string(seg.get_end(0).get_name()->get_str()) + ", ";
//    insert_str += base::quoted_string(seg.get_end_id(0)->get_str()) + ");";
//    conn.exec(insert_str);
//
//
//}

//void
//ResourceBuilder::_build_basic_library(
//        String const & path,
//        util::SegmentType seg_type,
//        std::map<String, int> const & excluded_motifs) {
//
//    auto lib_type = base::filename(path);
//    auto lib_path = base::resources_path() + "/motif_libraries/"+lib_type+".db";
//
//    LOGI << "building "+lib_type+" segement library at: " + lib_path;
//
//    try { std::remove(lib_path.c_str()); }
//    catch(...) { }
//
//    auto db = util::sqlite::Database(lib_path);
//    auto conn = util::sqlite::Connection(db);
//    util::sqlite::create_table(conn, motif_table_);
//
//    DIR *pDIR = opendir(path.c_str());
//    struct dirent *entry;
//    auto i = 0;
//
//    while ((entry = readdir(pDIR)) != NULL) {
//        auto fname = String(entry->d_name);
//        if(seg_type == util::SegmentType::HELIX) {
//            auto spl = base::split_str_by_delimiter(fname, ".");
//            if (spl[1] == "IDEAL" || spl[1] == "LE") { continue; }
//        }
//        //if(fname.length() < 10) { continue;}
//        if(fname[0] == '.') { continue; }
//        if(excluded_motifs.find(fname) != excluded_motifs.end()) { continue; }
//
//        auto pdb_path = path + fname + "/" + fname + ".pdb";
//        auto segs = seg_f_.all_segments_from_pdb(pdb_path, seg_type, false);
//
//        for(auto & seg : segs ) {
//            seg_f_.align_segment_to_ref_frame(*seg);
//            _insert_segment_to_motif_table(*seg, i, conn);
//            i++;
//        }
//
//    }
//
//    closedir(pDIR);
//    delete entry;
//}



}























