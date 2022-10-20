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

    _libnames = get_libnames();
    _rng = util::RandomNumberGenerator();
    auto path = _get_path(libname);
    SegmentSqliteConnection conn(path);
    _connection = conn;
    _max_size = _connection.count();
    // max_size_ = 1; // Does this matter? CJ
  }

  ~SegmentSqliteLibrary() {}

public: // iterator stuff
  class iterator {
  public:
    iterator(std::map<String, structure::SegmentOP>::iterator const &i)
        : _i(i) {}

    iterator operator++() {
      _i++;
      return *this;
    }

    structure::SegmentOP const &operator*() { return _i->second; }

    bool operator==(iterator const &rhs) const { return _i == rhs._i; }

    bool operator!=(iterator const &rhs) const { return _i != rhs._i; }

  private:
    std::map<String, structure::SegmentOP>::iterator _i;
  };

  iterator begin() { return iterator(_data.begin()); }

  iterator end() { return iterator(_data.end()); }

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
  SegmentSqliteConnection _connection;
  std::map<String, structure::SegmentOP> _data;
  util::RandomNumberGenerator _rng;
};

typedef std::shared_ptr<SegmentSqliteLibrary> SegmentSqliteLibraryOP;

} // namespace resources

#endif /* defined(__RNAMake__segment_sqlite_library__) */
