//
//  residue.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "secondary_structure/residue.h"

namespace secondary_structure {

ResType
convert_res_name_to_type(
        char c) {
    if (c == 'A')      { return ResType::A; }
    else if (c == 'C') { return ResType::C; }
    else if (c == 'G') { return ResType::G; }
    else if (c == 'U') { return ResType::U; }
    else if (c == 'T') { return ResType::U; }
    else if (c == 'W') { return ResType::W; }
    else if (c == 'S') { return ResType::S; }
    else if (c == 'M') { return ResType::M; }
    else if (c == 'K') { return ResType::K; }
    else if (c == 'R') { return ResType::R; }
    else if (c == 'Y') { return ResType::Y; }
    else if (c == 'B') { return ResType::B; }
    else if (c == 'D') { return ResType::D; }
    else if (c == 'H') { return ResType::H; }
    else if (c == 'V') { return ResType::V; }
    else if (c == 'N') { return ResType::N; }
    else {
        throw secondary_structure::Exception("incorrect character for secondary string");
    }
}

String
convert_res_type_to_str(
        ResType r) {
    if(r == ResType::A)      { return "A"; }
    else if(r == ResType::C) { return "C"; }
    else if(r == ResType::G) { return "G"; }
    else if(r == ResType::U) { return "U"; }
    else if(r == ResType::W) { return "W"; }
    else if(r == ResType::S) { return "S"; }
    else if(r == ResType::M) { return "M"; }
    else if(r == ResType::K) { return "K"; }
    else if(r == ResType::R) { return "R"; }
    else if(r == ResType::Y) { return "Y"; }
    else if(r == ResType::B) { return "B"; }
    else if(r == ResType::D) { return "D"; }
    else if(r == ResType::H) { return "H"; }
    else if(r == ResType::V) { return "V"; }
    else if(r == ResType::N) { return "N"; }
    else {
        throw secondary_structure::Exception("incorrect ResType cannot convert to string: " + std::to_string((int)r));
    }

}

bool
does_restype_satisfy_constraint(
        ResType r,
        ResType constraint) {
    if(constraint == ResType::A)      { return r == ResType::A; }
    else if(constraint == ResType::C) { return r == ResType::C; }
    else if(constraint == ResType::G) { return r == ResType::G; }
    else if(constraint == ResType::U) { return r == ResType::U; }
    else if(constraint == ResType::W) { return is_restype_a_weak(r); }
    else if(constraint == ResType::S) { return is_restype_a_strong(r); }
    else if(constraint == ResType::M) { return is_restype_a_amino(r); }
    else if(constraint == ResType::K) { return is_restype_a_keto(r); }
    else if(constraint == ResType::R) { return is_restype_a_purine(r); }
    else if(constraint == ResType::Y) { return is_restype_a_pyrimidine(r); }
    else if(constraint == ResType::B) { return is_restype_not_A(r); }
    else if(constraint == ResType::D) { return is_restype_not_C(r); }
    else if(constraint == ResType::H) { return is_restype_not_G(r); }
    else if(constraint == ResType::V) { return is_restype_not_U(r); }
    else if(constraint == ResType::N) { return true; }
    else {
        throw secondary_structure::Exception("unknown constraint ResType: " + std::to_string((int)constraint));
    }
}

bool
is_restype_a_weak(
        ResType r) {
    return (r == ResType::A || r == ResType::U);

}

bool
is_restype_a_strong(
        ResType r) {
    return (r == ResType::C || r == ResType::G);
}

bool
is_restype_a_amino(
        ResType r) {
    return (r == ResType::A || r == ResType::C);
}

bool
is_restype_a_keto(
        ResType r) {
    return (r == ResType::G || r == ResType::U);
}

bool
is_restype_a_purine(
        ResType r) {
    return (r == ResType::A || r == ResType::G);

}

bool
is_restype_a_pyrimidine(
        ResType r) {
    return (r == ResType::C || r == ResType::U);

}

bool
is_restype_not_A(
        ResType r) {
    return (r == ResType::C || r == ResType::G || r == ResType::U);
}

bool
is_restype_not_C(
        ResType r) {
    return (r == ResType::A || r == ResType::G || r == ResType::U);
}

bool
is_restype_not_G(
        ResType r) {
    return (r == ResType::A || r == ResType::C || r == ResType::U);
}

bool
is_restype_not_U(
        ResType r) {
    return (r == ResType::A || r == ResType::C || r == ResType::G);
}

bool
is_restype_a_ambiguous_code(
        ResType r) {
    return ! (r == ResType::A || r == ResType::C || r == ResType::G || r == ResType::U);
}


}
