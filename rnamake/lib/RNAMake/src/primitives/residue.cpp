//
// Created by Joseph Yesselman on 1/26/17.
//

#include "primitives/residue.h"

namespace primitives {

Residue::Residue(
        String const & name,
        int num,
        String const & chain_id):
        name_(name),
        num_(num),
        chain_id_(chain_id),
        i_code_(String("")),
        uuid_(Uuid()) {}

Residue::Residue(
        String const & name,
        int num,
        String const & chain_id,
        String const & i_code):
        name_(name),
        num_(num),
        chain_id_(chain_id),
        i_code_(i_code),
        uuid_(Uuid()) {}

Residue::Residue(
        String const & name,
        int num,
        String const & chain_id,
        String const & i_code,
        Uuid const & uuid):
        name_(name),
        num_(num),
        chain_id_(chain_id),
        i_code_(i_code),
        uuid_(uuid) {}
}
