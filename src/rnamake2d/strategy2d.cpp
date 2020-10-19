#include <rnamake2d/strategy2d.h>

#include <rnamake2d/strategy/bad_helix3.h>
#include <rnamake2d/strategy/bad_single_strand3.h>
#include <rnamake2d/strategy/bad_hairpin4_doublet.h>
#include <rnamake2d/strategy/pct_gc_pairs.h>
#include <rnamake2d/strategy/bad_helix2.h>
#include <rnamake2d/strategy/helix1.h>
#include <rnamake2d/strategy/junction_closing.h>

namespace rnamake2d {

    Strategy2DOP
    get_strategy( String const & strat ) {
        if(strat == "BadHelix3") {
            return std::make_shared<rnamake2d::BadHelix3>();
        } else if (strat == "BadSingleStrand3") {
            return std::make_shared<rnamake2d::BadSingleStrand3>();
        } else if (strat == "BadHairpin4Doublet") {
            return std::make_shared<rnamake2d::BadHairpin4Doublet>();
        } else if (strat == "PctGCPairs") {
            return std::make_shared<rnamake2d::PctGCPairs>();
        } else if (strat == "BadHelix2") {
            return std::make_shared<rnamake2d::BadHelix2>();
        } else if (strat == "Helix1") {
            return std::make_shared<rnamake2d::Helix1>();
        } else if (strat == "JunctionClosing") {
            return std::make_shared<rnamake2d::JunctionClosing>();
        }
    }
}