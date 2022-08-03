//
//  motif_state_ensemble_sqlite_library.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_ensemble_sqlite_library__
#define __RNAMake__motif_state_ensemble_sqlite_library__

#include <stdio.h>

#include <stdio.h>

#include "motif/motif_state_ensemble.h"
#include "resources/motif_ensemble_sqlite_connection.h"
#include "resources/motif_sqlite_library.h"
#include "util/random_number_generator.h"

namespace resources {

class MotifStateEnsembleSqliteLibrary : public SqliteLibrary {
public:
  MotifStateEnsembleSqliteLibrary(String const &libname) {

    _libnames = get_libnames();
    _rng = util::RandomNumberGenerator();
    auto path = _get_path(libname);
    MotifEnsembleSqliteConnection conn(path);
    _connection = conn;
    _max_size = _connection.count();
  }

  ~MotifStateEnsembleSqliteLibrary() {}

public: // iterator stuff
  class iterator {
  public:
    iterator(std::map<String, motif::MotifStateEnsembleOP>::iterator const &i)
        : _i(i) {}

    iterator operator++() {
      ++_i;
      return *this;
    }

    motif::MotifStateEnsembleOP const &operator*() { return _i->second; }

    bool operator==(iterator const &rhs) const { return _i == rhs._i; }

    bool operator!=(iterator const &rhs) const { return _i != rhs._i; }

  private:
    std::map<String, motif::MotifStateEnsembleOP>::iterator _i;
  };

  iterator begin() { return iterator(_data.begin()); }

  iterator end() { return iterator(_data.end()); }

public:
  static StringStringMap get_libnames();

  motif::MotifStateEnsembleOP get(String const &name = dummy_name,
                                  String const &id = dummy_id);

  motif::MotifStateEnsembleOPs get_multi(String const &name = dummy_name,
                                         String const &id = dummy_id);

  int contains(String const &name = dummy_name, String const &id = dummy_id);

  motif::MotifStateEnsembleOP get_random();

  void load_all(int limit = 99999);

private:
  String _generate_query(String const &, String const &);

private:
  MotifEnsembleSqliteConnection _connection;
  std::map<String, motif::MotifStateEnsembleOP> _data;
  util::RandomNumberGenerator _rng;
};

typedef std::shared_ptr<MotifStateEnsembleSqliteLibrary>
    MotifStateEnsembleSqliteLibraryOP;

} // namespace resources

#endif /* defined(__RNAMake__motif_state_ensemble_sqlite_library__) */
