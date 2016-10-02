
//headers for testing
#include "../common.hpp"
#include "../tools/motif_graph_builder.hpp"


//RNAMake Headers
#include "base/settings.h"
#include "math/numerical.h"
#include "motif/motif.h"
#include "motif/motif_factory.h"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_topology.h"
#include "sequence_optimizer/sequence_optimizer_3d.hpp"

TEST_CASE( "Test Sequence Optimizer", "[SequenceOptimizer]" ) {
    
    auto builder = MotifGraphBuilder();
    auto mg = builder.build(2);
    auto m = RM::instance().motif("HAIRPIN.1C0A.0");
    mg->add_motif(m);
    mg->replace_ideal_helices();
    
    auto c = GraphtoTree();
    auto mt = c.convert(mg);
    
    auto so = SequenceOptimizer3D();
    auto sols = so.get_optimized_sequences(mt, mt->last_node()->data()->ends()[0],
                                           mt->last_node()->index(), 0);
    
    
    
}