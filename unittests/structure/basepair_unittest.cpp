

#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "structure/structure.h"
#include "structure/basepair.h"
#include "structure/is_equal.h"


TEST_CASE( "Test Basepairs for Structure" ) {
    auto path = base::unittest_resource_dir() + "/structure/test_str_to_structure.dat";
    auto lines =base::get_lines_from_file(path);
    auto rts = structure::ResidueTypeSet();
    auto s = std::make_shared<structure::Structure>(lines[0], rts);
    
    SUBCASE("Test creation from resiudes") {
        auto res1 = s->get_residue(103, "A", "");
        auto res2 = s->get_residue(104, "A", "");
        auto r = math::Matrix(0.0);

        auto bp = std::make_shared<structure::Basepair>(res1, res2, r, "c...");
        
        CHECK(bp->res1() == res1);
        CHECK(bp->res2() == res2);
    }

}
