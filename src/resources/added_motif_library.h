//
//  added_motif_library.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/9/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__added_motif_library__
#define __RNAMake__added_motif_library__

#include <stdio.h>

#include "motif/motif.h"
#include "resources/motif_sqlite_library.h"

namespace resources {

/*
 * Exception for added motif library
 */
class AddedMotifLibraryException : public std::runtime_error {
public:
  /**
   * Standard constructor for AddedMotifLibraryException
   * @param   message   Error message for added motif library
   */
  AddedMotifLibraryException(String const &message)
      : std::runtime_error(message) {}
};

class AddedMotifLibrary {
public:
  AddedMotifLibrary() : motifs_(motif::MotifOPs()) {}

  ~AddedMotifLibrary() {}

public:
  void add_motif(motif::MotifOP const &m) {
    if (_find_motifs(m->name(), m->end_ids()[0], m->ends()[0]->name()).size() !=
        0) {
      throw AddedMotifLibraryException(
          "trying to add the same motif twice to library");
    }

    motifs_.push_back(m);
  }

  motif::MotifOP get(String const &name = dummy_name,
                     String const &end_id = dummy_end_id,
                     String const &end_name = dummy_name);

  motif::MotifOPs get_multi(String const &name = dummy_name,
                            String const &end_id = dummy_end_id,
                            String const &end_name = dummy_name);

  int contains(String const &name = dummy_name,
               String const &end_id = dummy_end_id,
               String const &end_name = dummy_name);

private:
  motif::MotifOPs _find_motifs(String const &name = dummy_name,
                               String const &end_id = dummy_end_id,
                               String const &end_name = dummy_name);

private:
  motif::MotifOPs _motifs;
};

} // namespace resources

#endif /* defined(__RNAMake__added_motif_library__) */
