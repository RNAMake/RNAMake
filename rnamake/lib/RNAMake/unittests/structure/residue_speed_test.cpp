
#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "math/numerical.h"
#include "util/x3dna.h"
#include "structure/residue.h"
//#include "structure/structure.h"
//#include "structure/basepair.h"
//#include "structure/is_equal.hpp"

state::Residue const &
copy_ref(state::ResidueOP const & rs) {
    return *rs;
}

TEST_CASE( "Test Speed of reidues", "[Residues]" ) {
    auto rts = ResidueTypeSet();
    auto path = unittest_resource_dir() + "residue/test_str_to_residue.dat";
    auto lines = get_lines_from_file(path);
    auto residues  = ResidueOPs();
    for(auto const & l : lines) {
        if(l.size() < 10) { break; } // end of file
        auto r = std::make_shared<Residue>(l, rts);
        residues.push_back(r);
    }

    /*for(int i = 0; i < 10000000; i++) {
        auto bp_copy = Basepair(*bp);
        bp_copy.move(Point(rand(), rand(), rand()));
    }*/
    residues[0]->build_beads();
    Matrices rots(100);
    Points ps(100);
    for(int i = 0; i < 100; i++) {
        rots[i] = Matrix(rand(), rand(), rand(),
                         rand(), rand(), rand(),
                         rand(), rand(), rand());
        ps[i] = Point(rand(), rand(), rand());
    }

    auto rs = residues[0]->get_state();
    auto dummy = Point();
    int pos = 0;
    for(int i = 0; i < 100000000; i++) {
        auto * r_copy = new state::Residue(*rs);
        pos = rand() % 100;
        r_copy->fast_transform(rots[pos], ps[pos], dummy);
        delete r_copy;

    }

}