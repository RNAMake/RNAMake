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
#include <rnamake2d/strategy/helix2.h>
#include <rnamake2d/strategy/non_canonical.h>
#include <rnamake2d/strategy/junction_non_gc_closing.h>
#include <rnamake2d/strategy/bad_juction3_3_3_3.h>
#include <rnamake2d/strategy/good_juction3_3_3_3.h>
#include <rnamake2d/strategy/good_tetraloop_doublet.h>
#include <rnamake2d/strategy/bad_junction2_1_1.h>
#include <rnamake2d/strategy/good_junction2_1_1.h>
#include <rnamake2d/strategy/bad_singlestrand4.h>
#include <rnamake2d/strategy/singlestrand5.h>
#include <rnamake2d/strategy/good_junction3_1_1_1.h>
#include <rnamake2d/strategy/good_tetraloop_triplet.h>

// eternabot ones
#include <rnamake2d/strategy/merryskies_only_as_in_the_loops.h>
#include <rnamake2d/strategy/merryskies_1_1_loop.h>
#include <rnamake2d/strategy/clollin_gs_in_place.h>
#include <rnamake2d/strategy/eli_no_blue_nucleotides_in_hook.h>
#include <rnamake2d/strategy/eli_double_AUPair_strategy.h>
#include <rnamake2d/strategy/penguian_clean_dotplot.h>
#include <rnamake2d/strategy/aldo_mismatch.h>
#include <rnamake2d/strategy/eli_red_line.h>
#include <rnamake2d/strategy/aldo_loops_and_stacks.h>

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
        } else if (strat == "Helix2") {
            return std::make_shared<rnamake2d::Helix2>();
        } else if (strat == "NonCanonical") {
            return std::make_shared<rnamake2d::NonCanonical>();
        } else if (strat == "JunctionNonGCClosing") {
            return std::make_shared<rnamake2d::JunctionNonGCClosing>();
        } else if (strat == "BadJunction3_3_3_3") {
            return std::make_shared<rnamake2d::BadJunction3_3_3_3>();
        } else if (strat == "GoodJunction3_3_3_3") {
            return std::make_shared<rnamake2d::GoodJunction3_3_3_3>();
        } else if (strat == "GoodTetraloopDoublet") {
            return std::make_shared<rnamake2d::GoodTetraloopDoublet>();
        } else if (strat == "BadJunction2_1_1") {
            return std::make_shared<rnamake2d::BadJunction2_1_1>();
        } else if (strat == "GoodJunction2_1_1") {
            return std::make_shared<rnamake2d::GoodJunction2_1_1>();
        } else if (strat == "BadSingleStrand4") {
            return std::make_shared<rnamake2d::BadSingleStrand4>();
        } else if (strat == "SingleStrand5") {
            return std::make_shared<rnamake2d::SingleStrand5>();
        } else if (strat == "GoodJunction3_1_1_1") {
            return std::make_shared<rnamake2d::GoodJunction3_1_1_1>();
        } else if (strat == "GoodTetraloopTriplet") {
            return std::make_shared<rnamake2d::GoodTetraloopTriplet>();
        } else if (strat == "merryskies_only_as_in_the_loops") {
            return std::make_shared<rnamake2d::MerrySkiesOnlyAsInTheLoops>();
        } else if (strat == "clollin_gs_in_place") {
            return std::make_shared<rnamake2d::ClollinGsInPlace>();
        } else if (strat == "eli_no_blue_nucleotides_in_hook") {
            return std::make_shared<rnamake2d::EliNoBlueNucleotidesInHook>();
        } else if (strat == "merryskies_1_1_loop") {
            return std::make_shared<rnamake2d::MerrySkies_1_1_Loop>();
        } else if (strat == "eli_double_AUPair_strategy") {
            return std::make_shared<rnamake2d::EliDoubleAUPair>();
        } else if (strat == "penguian_clean_dotplot") {
            return std::make_shared<rnamake2d::PenguianCleanDotplot>();
        } else if (strat == "aldo_mismatch") {
            return std::make_shared<rnamake2d::AldoMismatch>();
        } else if (strat == "eli_red_line") {
            return std::make_shared<rnamake2d::EliRedLine>();
        } else if (strat == "aldo_loops_and_stacks") {
            return std::make_shared<rnamake2d::AldoLoopsAndStacks>();
        }
        else {
            throw base::RNAMakeException("ERROR: the strategy " + strat + " is not implemented yet.");
        }
    }
}