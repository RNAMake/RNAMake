#ifndef __RNAMAKE_SCORE_FUNCTION_H__
#define __RNAMAKE_SCORE_FUNCTION_H__

#include <map>

#include <base/types.h>
#include <rnamake2d/design.h>
#include <rnamake2d/Rule.h>
#include <eternabot/strategy.h>
#include <vienna/vienna.h>

#include <rnamake2d/strategy2d.h>

namespace rnamake2d {

    class ScoreFunction {
        eternabot::StrategyOPs strategies_;
        Reals weights_;
        int num_strats;
    public:
        ScoreFunction() {

            strategies_ = {

                get_strategy("BadHelix3"),
                get_strategy("BadSingleStrand3"),
                get_strategy("BadHairpin4Doublet"),
                get_strategy("PctGCPairs"),
                get_strategy("BadHelix2"),
                get_strategy("Helix1"),
                get_strategy("JunctionClosing"),
                get_strategy("SingleStrandTripleA"),
                get_strategy("BadHairpin4Triplet"),
                get_strategy("GRepeats"),
                get_strategy("GoodTetraloopSinglet"),
                get_strategy("Helix4"),
                get_strategy("BadSingleStrand5"),
                get_strategy("Helix3"),
                get_strategy("BadHelix4"),
                get_strategy("Helix2"),
                get_strategy("NonCanonical"),
                get_strategy("JunctionNonGCClosing"),
                get_strategy("BadJunction3_3_3_3"),
                get_strategy("GoodJunction3_3_3_3"),
                get_strategy("GoodTetraloopDoublet"),
                get_strategy("BadJunction2_1_1"),
                get_strategy("GoodJunction2_1_1"),
                get_strategy("GoodHairpin4"),
                get_strategy("BadSingleStrand4"),
                get_strategy("SingleStrand5"),
                get_strategy("GoodJunction3_1_1_1"),
                get_strategy("GoodTetraloopTriplet"),
                get_strategy("merryskies_only_as_in_the_loops"),
                get_strategy("clollin_gs_in_place"),
                get_strategy("eli_no_blue_nucleotides_in_hook"),
                get_strategy("merryskies_1_1_loop"),
                get_strategy("eli_double_AUPair_strategy"),
                get_strategy("penguian_clean_dotplot"),
                get_strategy("aldo_mismatch"),
                get_strategy("eli_red_line"),
                get_strategy("aldo_loops_and_stacks"),
                get_strategy("eli_green_blue_strong_middle_half"),
                get_strategy("deivad_deivad_strategy"),
                get_strategy("dejerpha_basic_test"),
                get_strategy("eli_twisted_basepairs"),
                get_strategy("eli_green_line"),
                get_strategy("example_gc60"),
                get_strategy("eli_direction_of_gc_pairs_in_multiloops_neckarea"),
                get_strategy("eli_multiloop_similarity"),
                get_strategy("ding_quad_energy"),
                get_strategy("berex_berex_loop_basic"),
                get_strategy("eli_legal_placement_of_GUpairs"),
                get_strategy("eli_wrong_direction_of_gc_pairs_in_multiloops"),
                get_strategy("eli_direction_of_gc_pairs_in_multiloops"),
                get_strategy("eli_no_blue_nucleotides_strategy"),
                get_strategy("berex_basic_test"),
                get_strategy("kkohli_test_by_kkohli"),
                get_strategy("AldoRepetition"),
                get_strategy("CJGreenLine"),
                get_strategy("CJMismatch"),

            } ;
            // weights
            weights_ = {

                    0.10845716043121478,    // BadHelix3
                    0.1902229692459677,     // BadSingleStrand3
                    0.11621958146693934,    // BadHairpin4Doublet
                    0.6973340034106996,     // PctGCPairs
                    0.051953978499975795,   // BadHelix2
                    0.5804069633989485,     // Helix1
                    -0.7514795405184205,    // JunctionClosing
                    -0.2524611477466161f,   // SingleStrandTripleA
                    0.14470067211618373f,   // BadHairpin4Triplet
                    -0.05032464330299019,   // GRepeats
                    0.8372216739418228,     // GoodTetraloopSinglet
                    0.6596505460982536,     // Helix4
                    0.14447503330605593,    // BadSingleStrand5
                    0.29379550128941323,    // Helix3
                    0.17028695543216954,    // BadHelix4
                    0.331594773993627,      // Helix2
                    -0.0853402523094039,    // NonCanonical
                    0.30753499424534314,    // JunctionNonGCClosing
                    -0.097400093245612,     // BadJunction3_3_3_3
                    1.1734967921515913,     // GoodJunction3_3_3_3
                    1.393824377152868,      // GoodTetraloopDoublet
                    0.24122961401641185,    // BadJunction2_1_1
                    -0.4598285526297501,    // GoodJunction2_1_1
                    -1.0079078237755852,    // GoodHairpin4
                    0.05046284982294122,    // BadSingleStrand4
                    -0.11077156860161347,   // SingleStrand5
                    -9.462788548128273,     // GoodJunction3_1_1_1
                    1.916441326283273,      // GoodTetraloopTriplet
                    -0.10825825690698382,   // merryskies_only_as_in_the_loops
                    0.2800622250719158,     // clollin_gs_in_place
                    0.1046622560893758,     // eli_no_blue_nucleotides_in_hook
                    0.16505853542440824,    // merryskies_1_1_loop
                    -0.19378279913278135,   // eli_double_AUPair_strategy
                    -0.8106666177952848,    // penguian_clean_dotplot
                    -0.19581776453734906,   // aldo_mismatch
                    -0.0018275413712541821, // eli_red_line
                    0.15067172153854896,    // aldo_loops_and_stacks,
                    -0.06418945089333378,   // eli_green_blue_strong_middle_half
                    0.08548851359538262,    // deivad_deivad_strategy
                    2.8236531196192134,     // dejerpha_basic_test
                    0.21601775412216717,    // eli_twisted_basepairs
                    0.007680977135075769,   // eli_green_line
                    -0.1489280364856781,    // example_gc60
                    0.09865127259968665,    // eli_direction_of_gc_pairs_in_multiloops_neckarea
                    0.00978877166914174,    // eli_multiloop_similarity
                    0.00978877166906401,    // ding_quad_energy
                    0.019145505229272428,   // berex_berex_loop_basic
                    0.05120674277182907,    // eli_legal_placement_of_GUpairs
                    -0.19731059143533014,   // eli_wrong_direction_of_gc_pairs_in_multiloops
                    0.09865127290814124,    // eli_direction_of_gc_pairs_in_multiloops
                    -0.8253549680053338,    // eli_no_blue_nucleotides_strategy
                    -0.2751340813578167,    // berex_basic_test
                    -0.06803891875360141,   // kkohli_test_by_kkohli
                    0.15268020542080335,    // AldoRepetition
                    -0.004643042781698109,  // EliGreenLine
                    0.2430050990096834,     // CJMismatch
            };


            num_strats = strategies_.size();
        }

        public:
            double
            score( Feature2DOP const & ) const ;


        private:
            std::map<double,RuleOP> scoring_rules_;
    };

    class ViennaScore {
        private:
            vienna::Vienna vienna_;
        public:
            double
            score(Design const&);

            double
            score_mutation(Design const&);
    };
}



#endif // __RNAMAKE_SCORE_FUNCTION_H__
// will hold off on these for now...
//NupackStructureComp,0.7965075352904527
//SpotStructureComp,0.4615538691139922

// => AldoRepetition,0.15268020542080335
// => EliGreenLine,-0.004643042781698109
//ViennaStructureComp,1.5272406678598347
//ViennaMFENormalized,-0.49663024818195145

// => CJMismatch,0.2430050990096834
//RnastructureStructureComp,-0.23831987534823496
// => merryskies_only_as_in_the_loops,-0.10825825690698382
// => dejerpha_basic_test,2.8236531196192134
// => clollin_gs_in_place,0.2800622250719158
// => eli_no_blue_nucleotides_in_hook,0.1046622560893758
// => merryskies_1_1_loop,0.16505853542440824
// => eli_double_AUPair_strategy,-0.19378279913278135
// => eli_green_blue_strong_middle_half,-0.06418945089333378
// => example_gc60,-0.1489280364856781
// => penguian_clean_dotplot,-0.8106666177952848
// => eli_twisted_basepairs,0.21601775412216717
// => aldo_loops_and_stacks,0.15067172153854896
// => eli_direction_of_gc_pairs_in_multiloops_neckarea,0.09865127259968665
// => eli_multiloop_similarity,0.00978877166914174
// => eli_green_line,0.007680977135075769
// => ding_quad_energy,0.00978877166906401
// => berex_berex_loop_basic,0.019145505229272428
// => eli_legal_placement_of_GUpairs,0.05120674277182907
// => aldo_mismatch,-0.19581776453734906
// => eli_red_line,-0.0018275413712541821
// => eli_wrong_direction_of_gc_pairs_in_multiloops,-0.19731059143533014
// => deivad_deivad_strategy,0.08548851359538262
// => eli_direction_of_gc_pairs_in_multiloops,0.09865127290814124
// => eli_no_blue_nucleotides_strategy,-0.8253549680053338
// => berex_basic_test,-0.2751340813578167
// => kkohli_test_by_kkohli,-0.06803891875360141
// => BadHelix3,0.10845716043121478
// => SingleStrandTripleA,-0.2524611477466161
// => BadHairpin4Doublet,0.11621958146693934
// => BadHelix2,0.051953978499975795
// => BadHairpin4Triplet,0.14470067211618373
// => BadSingleStrand3,0.1902229692459677
// => PctGCPairs,0.6973340034106996
// => Helix1,0.5804069633989485
// => JunctionClosing,-0.7514795405184205
// => SingleStrand3,0.08819793678377201
// => GRepeats,-0.05032464330299019
// => GoodTetraloopSinglet,0.8372216739418228
// => Helix4,0.6596505460982536
// => BadJunction2_1_1,0.24122961401641185
// => BadSingleStrand5,0.14447503330605593
// => Helix3,0.29379550128941323
// => GoodTetraloopDoublet,1.393824377152868
// => GoodJunction2_1_1,-0.4598285526297501
// => GoodJunction3_1_1_1,-9.462788548128273
// => GoodTetraloopTriplet,1.916441326283273
// => SingleStrand5,-0.11077156860161347
// => GoodJunction3_3_3_3,1.1734967921515913
// => Badjunction3_3_3_3,-0.097400093245612
// => BadSingleStrand4,0.05046284982294122
// => Helix2,0.331594773993627
// => GoodHairpin4,-1.0079078237755852
// => NonCanonical,-0.0853402523094039
// => BadHelix4,0.17028695543216954
// => JunctionNonGCClosing,0.30753499424534314