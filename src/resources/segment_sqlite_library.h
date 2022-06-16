//
//  segment_sqlite_library.h
//  RNAMake
//
//  Created by Hassan Abdelsamad on 8/26/21.
//

#ifndef __RNAMake__segment_sqlite_library__
#define __RNAMake__segment_sqlite_library__

#include <iostream>
#include <stdio.h>

#include "base/settings.h"
#include "resources/segment_sqlite_connection.h"
#include "resources/sqlite_library.h"
#include "structure/segment.h"
#include "util/random_number_generator.h"

namespace resources {

static String dummy_name = "", dummy_end_id = "", dummy_end_name = "",
              dummy_id = "";

class SegmentSqliteLibrary : public SqliteLibrary {
public:
  SegmentSqliteLibrary(String const &libname) {

    libnames_ = get_libnames();
    rng_ = util::RandomNumberGenerator();
    auto path = _get_path(libname);
    SegmentSqliteConnection conn(path);
    connection_ = conn;
    max_size_ = connection_.count();
    // max_size_ = 1; // Does this matter? CJ
  }

  ~SegmentSqliteLibrary() {}

public: // iterator stuff
  class iterator {
  public:
    iterator(std::map<String, structure::SegmentOP>::iterator const &i)
        : i_(i) {}

    iterator operator++() {
      i_++;
      return *this;
    }

    structure::SegmentOP const &operator*() { return i_->second; }

    bool operator==(iterator const &rhs) const { return i_ == rhs.i_; }

    bool operator!=(iterator const &rhs) const { return i_ != rhs.i_; }

  private:
    std::map<String, structure::SegmentOP>::iterator i_;
  };

  iterator begin() { return iterator(data_.begin()); }

  iterator end() { return iterator(data_.end()); }

public:
  static StringStringMap get_libnames();

  structure::SegmentOP get(String const &name = dummy_name,
                           String const &end_id = dummy_end_id,
                           String const &end_name = dummy_name,
                           String const &id = dummy_id);

  structure::SegmentOPs get_multi(String const &name = dummy_name,
                                  String const &end_id = dummy_end_id,
                                  String const &end_name = dummy_name,
                                  String const &id = dummy_id);

  int contains(String const &name = dummy_name,
               String const &end_id = dummy_end_id,
               String const &end_name = dummy_name,
               String const &id = dummy_id);

  structure::SegmentOP get_random();

  void load_all(int limit = 99999);

private:
  String _generate_query(String const &, String const &, String const &,
                         String const &);

private:
  SegmentSqliteConnection connection_;
  std::map<String, structure::SegmentOP> data_;
  util::RandomNumberGenerator rng_;
};

typedef std::shared_ptr<SegmentSqliteLibrary> SegmentSqliteLibraryOP;

} // namespace resources

#endif /* defined(__RNAMake__segment_sqlite_library__) */
