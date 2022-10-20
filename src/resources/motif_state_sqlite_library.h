//
//  motif_state_sqlite_library.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_sqlite_library__
#define __RNAMake__motif_state_sqlite_library__

#include <stdio.h>

//#include "motif/motif_state.h"
#include "resources/motif_sqlite_connection.h"
#include "resources/motif_sqlite_library.h"
#include "util/random_number_generator.h"

namespace resources {

class MotifStateSqliteLibrary : public SqliteLibrary {
public:
  MotifStateSqliteLibrary(String const &libname) {

    _libnames = get_libnames();
    _rng = util::RandomNumberGenerator();
    _name = libname;
    auto path = _get_path(libname);
    MotifSqliteConnection conn(path);
    _connection = conn;
    _max_size = _connection.count();
  }

  ~MotifStateSqliteLibrary() {}

public: // iterator stuff
  class iterator {
  public:
    iterator(std::map<String, motif::MotifStateOP>::iterator const &i)
        : _i(i) {}

    iterator operator++() {
      _i++;
      return *this;
    }

    motif::MotifStateOP const &operator*() { return _i->second; }

    bool operator==(iterator const &rhs) const { return _i == rhs._i; }

    bool operator!=(iterator const &rhs) const { return _i != rhs._i; }

  private:
    std::map<String, motif::MotifStateOP>::iterator _i;
  };

  iterator begin() { return iterator(_data.begin()); }

  iterator end() { return iterator(_data.end()); }

public:
  static StringStringMap get_libnames();

  motif::MotifStateOP get(String const &name = dummy_name,
                          String const &end_id = dummy_end_id,
                          String const &end_name = dummy_name,
                          String const &id = dummy_id);

  motif::MotifStateOPs get_multi(String const &name = dummy_name,
                                 String const &end_id = dummy_end_id,
                                 String const &end_name = dummy_name,
                                 String const &id = dummy_id);

  int contains(String const &name = dummy_name,
               String const &end_id = dummy_end_id,
               String const &end_name = dummy_name,
               String const &id = dummy_id);

  motif::MotifStateOP get_random();

  void load_all(int limit = 99999);

  String const &get_name() { return _name; }

private:
  String _generate_query(String const &, String const &, String const &,
                         String const &);

private:
  MotifSqliteConnection _connection;
  std::map<String, motif::MotifStateOP> _data;
  util::RandomNumberGenerator _rng;
  String _name;
};

typedef std::shared_ptr<MotifStateSqliteLibrary> MotifStateSqliteLibraryOP;

} // namespace resources

#endif /* defined(__RNAMake__motif_state_sqlite_library__) */
