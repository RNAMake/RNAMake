
#include <util/random_number_generator.h>
#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "math/numerical.h"
#include "util/x3dna.h"
#include "structure/structure.h"
#include "structure/basepair.h"

state::BasepairOP
_get_bp_from_str(String const & s) {
    auto spl = split_str_by_delimiter(s, ";");
    auto d = vector_from_str(spl[0]);
    auto r = matrix_from_str(spl[1]);
    auto sugars = vectors_from_str(spl[2]);
    auto bp_name_str = String("test");
    auto bp_name = std::make_shared<SimpleString>(bp_name_str);
    return std::make_shared<state::Basepair>(Uuid(), Uuid(), r, d, sugars[0], sugars[1], bp_name,
                                             X3dna::X3dnaBPType::cDDD,
                                             primitives::Basepair::BasepairType::WC, Uuid());
}

TEST_CASE( "Test Speed of basepairs", "[Basepair]" ) {
    auto bps = state::BasepairOPs();
    auto path = unittest_resource_dir() + "/structure/test_get_transformed_state.dat";
    auto lines = get_lines_from_file(path);
    for(auto const & l : lines) {
        if (l.length() < 5) { break; }
        auto spl = split_str_by_delimiter(l, "|");
        auto bp1 = _get_bp_from_str(spl[0]);
        auto bp2 = _get_bp_from_str(spl[1]);
        bps.push_back(bp1);
        bps.push_back(bp2);
    }

    auto transform_info = state::TransformInfo();
    auto rng = RandomNumberGenerator();

    auto bp1 = state::BasepairOP(nullptr);
    auto bp2 = state::BasepairOP(nullptr);

    for(int i = 0; i < 10000000; i++) {
        bp1 = std::make_shared<state::Basepair>(*bps[rng.randrange(bps.size()-1)]);
        bp2 = std::make_shared<state::Basepair>(*bps[rng.randrange(bps.size()-1)]);

        bp1->get_transforming_r_and_t(*bp2, transform_info);
        bp2->fast_transform(transform_info);
    }


}