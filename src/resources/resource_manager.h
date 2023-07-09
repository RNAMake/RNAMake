//
//  resource_manager.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__resource_manager__
#define __RNAMake__resource_manager__

#include <stdio.h>

//RNAMake Headers
#include "base/types.h"
#include "motif/motif_factory.h"
#include "motif/motif_state.h"
#include "motif/motif_ensemble.h"
#include "resources/motif_sqlite_library.h"
#include "resources/motif_state_sqlite_library.h"
#include "resources/motif_state_ensemble_sqlite_library.h"
#include "resources/added_motif_library.h"

namespace resources {

class ResourceManagerException : public std::runtime_error {
public:
  ResourceManagerException(
    String const & message) :
    std::runtime_error(message) {
  }
};

class Manager { //RM for ResourceManager
protected:

  Manager();

  Manager(Manager const &); //Prevent construction
  void
  operator =(Manager const &);

private:

  ~Manager() = default;

public:

  static Manager &
  instance() {
    static Manager instance;
    return instance;
  }

public: // getting functions
  motif::MotifOP
  bp_step(
    String const &);

  motif::MotifStateOP
  bp_step_state(
    String const &);

  motif::MotifOP
  motif(
    String const & name = dummy_name,
    String const & end_id = dummy_end_id,
    String const & end_name = dummy_name);

  motif::MotifStateOP
  motif_state(
    String const & name = dummy_name,
    String const & end_id = dummy_end_id,
    String const & end_name = dummy_name);

  motif::MotifStateEnsembleOP
  motif_state_ensemble(
    String const & name = dummy_name);

public: // adding functions

  structure::RNAStructureOP
  get_structure(
    String const & path,
    String name = "",
    int force_num_chains = -1);

  void
  add_motif(
    String const & path,
    String name = "",
    util::MotifType mtype = util::MotifType::UNKNOWN);

  void
  add_motif(
    motif::MotifOP const & m,
    String name = "");

  void
  register_motif(
    motif::MotifOP const &);

  void
  register_extra_motif_ensembles(
    String const &);

  void
  register_motif_ensemble(
    String const &,
    String const &,
    motif::MotifEnsembleOP const &);

  int
  has_supplied_motif_ensemble(
    String const &,
    String const &);

  std::vector<motif::MotifEnsembleOP>
  get_supplied_motif_ensembles(
    String const &);

  motif::MotifEnsembleOP const &
  get_supplied_motif_ensemble(
    String const &,
    String const &);

public: // new functions that I would like in new code

  motif::MotifOP
  get_motif_from_state(
    motif::MotifStateOP ms) {
    auto m = motif(ms->name(), "", ms->end_names()[0]);
    align_motif(ms->end_states()[0], m->ends()[0], m);
    return m;
  }

private:
  std::map<String, MotifSqliteLibraryOP> mlibs_;
  std::map<String, MotifStateSqliteLibraryOP> ms_libs_;
  std::map<String, MotifStateEnsembleSqliteLibraryOP> mse_libs_;
  std::map<String, motif::MotifEnsembleOP> extra_me_;
  motif::MotifFactory mf_;
  AddedMotifLibrary added_motifs_;

};

inline
motif::MotifOP
get_motif_from_resource_manager(
  String const & name = dummy_name,
  String const & end_id = dummy_end_id,
  String const & end_name = dummy_name) {

  return Manager::instance().motif(name, end_id, end_name);

}

}

#endif /* defined(__RNAMake__resource_manager__) */
