//
// Created by Joseph Yesselman on 1/28/17.
//

#include "chain.h"

namespace state {

String
Chain::to_str() {
    auto s =  String("");
    for(auto const & r : residues_) {
        s += r->to_str() + ";";
    }
    return s;
}

}