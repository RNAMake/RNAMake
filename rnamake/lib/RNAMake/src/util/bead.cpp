//
// Created by Joseph Yesselman on 1/26/17.
//

#include "util/bead.h"

String
Bead::to_str() const{
    return vector_to_str(center_) + "," + std::to_string(btype_);
}
