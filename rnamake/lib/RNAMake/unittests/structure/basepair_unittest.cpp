
#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "structure/structure.h"
#include "structure/basepair.h"
#include "structure/is_equal.hpp"


TEST_CASE( "Test Basepairs for Structure", "[Basepair]" ) {
    auto path = unittest_resource_dir() + "/structure/test_str_to_structure.dat";
    auto lines =base::get_lines_from_file(path);
    auto rts = ResidueTypeSet();
    auto s = std::make_shared<Structure>(lines[0], rts);
    
    SECTION("Test creation from resiudes") {
        auto res1 = s->get_residue(103, "A", "");
        auto res2 = s->get_residue(104, "A", "");
        auto r = Matrix(0.0);

        auto bp = std::make_shared<Basepair>(res1, res2, r, "c...");
        
        REQUIRE(bp->res1() == res1);
        REQUIRE(bp->res2() == res2);
    }

}