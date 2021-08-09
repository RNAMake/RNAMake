//
// Created by Joseph Yesselman on 1/5/18.
//

#ifndef RNAMAKE_NEW_MOTIF_SQLITE_LIBRARY_H
#define RNAMAKE_NEW_MOTIF_SQLITE_LIBRARY_H

#include <structure/segment.h>
#include <resources/sqlite_library.h>

namespace resources {

class SegmentSqliteLibrary : public SqliteLibrary {
public:
    SegmentSqliteLibrary(
            String const & db_path,
            String const & table_name,
            structure::ResidueTypeSet const & rts):
            SqliteLibrary(db_path, table_name),
            rts_(rts),
            retrieved_columns_(Strings{"id", "data"}),
            segments_(std::map<int, structure::SegmentOP>()) {}

    ~SegmentSqliteLibrary() {}

public:
    structure::SegmentOP
    get_segment(
            StringStringMap const &) const;

    bool
    contains_segment(
            StringStringMap const &) const;

protected:
    Strings retrieved_columns_;
    structure::ResidueTypeSet const & rts_;
    mutable std::map<int, structure::SegmentOP> segments_; // ack this seems bad

};

}

#endif //RNAMAKE_NEW_MOTIF_SQLITE_LIBRARY_H
