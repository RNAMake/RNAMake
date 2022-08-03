//
//  motif_sqlite_library.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_sqlite_library__
#define __RNAMake__motif_sqlite_library__

#include <filesystem>
#include <iostream>
#include <stdio.h>

#include "base/settings.h"
#include "motif/motif.h"
#include "resources/motif_sqlite_connection.h"
#include "resources/sqlite_library.h"
#include "util/random_number_generator.h"

namespace resources {

static String dummy_name = "", dummy_end_id = "", dummy_end_name = "",
              dummy_id = "";

class MotifSqliteLibrary : public SqliteLibrary {
public:
  MotifSqliteLibrary(String const &libname) {

    _libnames = get_libnames();
    _rng = util::RandomNumberGenerator();
    auto path = _get_path(libname);
    MotifSqliteConnection conn(path);
    _connection = conn;
    _max_size = _connection.count();
    // max_size_ = 1; // Does this matter? CJ
  }
  // added by CJ for validation purposes
  MotifSqliteLibrary(int i, // something to change the overloading
                     std::filesystem::path const &path) {

    _libnames = get_libnames();
    _rng = util::RandomNumberGenerator();

    MotifSqliteConnection conn(path);
    _connection = conn;
    _max_size = _connection.count();
    // max_size_ = 1; // Does this matter? CJ
  }
  ~MotifSqliteLibrary() {}

public: // iterator stuff
  class iterator {
  public:
    iterator(std::map<String, motif::MotifOP>::iterator const &i) : _i(i) {}

    iterator operator++() {
      _i++;
      return *this;
    }

    motif::MotifOP const &operator*() { return _i->second; }

    bool operator==(iterator const &rhs) const { return _i == rhs._i; }

    bool operator!=(iterator const &rhs) const { return _i != rhs._i; }

  private:
    std::map<String, motif::MotifOP>::iterator _i;
  };

  iterator begin() { return iterator(_data.begin()); }

  iterator end() { return iterator(_data.end()); }

public:
  static StringStringMap get_libnames();

  motif::MotifOP get(String const &name = dummy_name,
                     String const &end_id = dummy_end_id,
                     String const &end_name = dummy_name,
                     String const &id = dummy_id);

  motif::MotifOPs get_multi(String const &name = dummy_name,
                            String const &end_id = dummy_end_id,
                            String const &end_name = dummy_name,
                            String const &id = dummy_id);

  int contains(String const &name = dummy_name,
               String const &end_id = dummy_end_id,
               String const &end_name = dummy_name,
               String const &id = dummy_id);

  motif::MotifOP get_random();

  void load_all(int limit = 99999);

private:
  static String _generate_query(String const &, String const &, String const &,
                                String const &);

private:
  MotifSqliteConnection _connection;
  std::map<String, motif::MotifOP> _data;
  util::RandomNumberGenerator _rng;
};

typedef std::shared_ptr<MotifSqliteLibrary> MotifSqliteLibraryOP;

} // namespace resources

#endif /* defined(__RNAMake__motif_sqlite_library__) */
