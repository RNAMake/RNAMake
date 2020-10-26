#include <rnamake2d/strategy2d.h>

#include <rnamake2d/strategy/bad_helix3.h>
#include <rnamake2d/strategy/bad_single_strand3.h>
#include <rnamake2d/strategy/bad_hairpin4_doublet.h>
#include <rnamake2d/strategy/pct_gc_pairs.h>
#include <rnamake2d/strategy/bad_helix2.h>
#include <rnamake2d/strategy/helix1.h>
#include <rnamake2d/strategy/junction_closing.h>
#include <rnamake2d/strategy/bad_0_1_bulge.h>
#include <rnamake2d/strategy/bad_hairpin4_triplet.h>
#include <rnamake2d/strategy/g_repeats.h>
#include <rnamake2d/strategy/good_hairpin4_singlet.h>
#include <rnamake2d/strategy/helix4.h>
#include <rnamake2d/strategy/bad_singlestrand5.h>
#include <rnamake2d/strategy/helix3.h>
#include <rnamake2d/strategy/bad_helix4.h>

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
        } else if (strat == "BadBulge_0_1") {
            return std::make_shared<rnamake2d::Bad_0_1_Bulge>();
//        } else if (strat == "BadSingleStrand3") {
//            return std::make_shared<rnamake2d::BadSingleStrand3>();
        } else if (strat == "SingleStrandTripleA") {
            return std::make_shared<rnamake2d::BadSingleStrand3>();
        } else if( strat == "BadHairpin4Triplet") {
            return std::make_shared<rnamake2d::BadHaripin4Triplet>();
        } else if (strat == "GRepeats") {
            return std::make_shared<rnamake2d::GRepeats>();
        } else if (strat == "GoodTetraloopSinglet") {
            return std::make_shared<rnamake2d::GoodHairpin4Singlet>();
        } else if (strat == "Helix4") {
            return std::make_shared<rnamake2d::Helix4>();
        } else if (strat == "BadSingleStrand5") {
            return std::make_shared<rnamake2d::BadSingleStrand5>();
        } else if (strat == "Helix3") {
            return std::make_shared<rnamake2d::Helix3>();
        } else if (strat == "BadHelix4") {
            return std::make_shared<rnamake2d::BadHelix4>();
        }
        else {
            throw base::RNAMakeException("ERROR: the strategy " + strat + " is not implemented yet.");
        }
    }
}