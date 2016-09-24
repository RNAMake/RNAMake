
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "math/numerical.h"
#include "motif/motif.h"
#include "motif/motif_factory.h"
#include "motif_tools/segmenter.h"

TEST_CASE( "Test Segmenting RNA Structures", "[Segmenter]" ) {
    auto mf = MotifFactory();
    auto path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6/";
    auto m = mf.motif_from_file(path);
    auto end1 = m->get_basepair("A129-A193")[0];
    auto end2 = m->get_basepair("A131-A192")[0];
    auto bps = BasepairOPs{end1, end2};
    
    auto segmenter = Segmenter();
    auto segments = segmenter.apply(m, bps);
    
    
    
    
}