#ifndef __RNAMAKE_SCORE_FUNCTION_H__
#define __RNAMAKE_SCORE_FUNCTION_H__

#include <map>

#include <base/types.h>
#include <rnamake2d/Design.h>
#include <rnamake2d/Rule.h>
#include <eternabot/strategy.h>
#include <vienna/vienna.h>

// actual strategies...

// first those from eternabot
#include <eternabot/strategy/berex_test.h>
#include <eternabot/strategy/clear_plot.h>
#include <eternabot/strategy/direction_of_gc.h>
#include <eternabot/strategy/modified_a_basic_test.h>
#include <eternabot/strategy/modified_berex_test.h>
#include <eternabot/strategy/modified_clear_plot.h>
#include <eternabot/strategy/modified_direction_of_gc.h>
#include <eternabot/strategy/modified_num_of_yellow.h>
#include <eternabot/strategy/num_of_yellow.h>

#include <rnamake2d/strategy2d.h>

namespace rnamake2d {

    class ScoreFunction {
        eternabot::StrategyOPs strategies_;
        Reals weights_;
        int num_strats;
    public:
        ScoreFunction() {

            strategies_ = {
                //std::make_shared<eternabot::BerexTest>(),
                std::make_shared<eternabot::CleanPlotStackCapsandSafeGC>(),
                std::make_shared<eternabot::DirectionofGCPairsinMultiLoops>(),
                std::make_shared<eternabot::ModifiedABasicTest>(),
                std::make_shared<eternabot::ModifiedBerexTest>(),
                //std::make_shared<eternabot::ModifiedNumofYellowNucleotidesperLengthofString>(),
                std::make_shared<eternabot::NumofYellowNucleotidesperLengthofString>(),
                get_strategy("BadHelix3"),
                get_strategy("BadSingleStrand3"),
                get_strategy("BadHairpin4Doublet"),
                get_strategy("PctGCPairs"),
                get_strategy("BadHelix2"),
                get_strategy("Helix1"),
                get_strategy("JunctionClosing")
            } ;
            // weights
            weights_ = {
                    0.1250677,
                    0.2156337,
                    0.09281782,
                    0.3661276,
                    0.2230357,
                    0.05,
                    0.10845716043121478,    // BadHelix3
                    0.1902229692459677,     // BadSingleStrand3
                    0.11621958146693934,    // BadHairpin4Doublet
                    0.6973340034106996,     // PctGCPairs
                    0.051953978499975795,   // BadHelix2
                    0.5804069633989485,     // Helix1
                    -0.7514795405184205     // JunctionClosing
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



#endif // __RNAMAKE_SCORE_FUNCTION_H__src/eternabot/strategy/a_basic_test.h


//merryskies_only_as_in_the_loops,-0.10825825690698382
//dejerpha_basic_test,2.8236531196192134
//clollin_gs_in_place,0.2800622250719158
//eli_no_blue_nucleotides_in_hook,0.1046622560893758
//merryskies_1_1_loop,0.16505853542440824
//eli_double_AUPair_strategy,-0.19378279913278135
//eli_green_blue_strong_middle_half,-0.06418945089333378
//example_gc60,-0.1489280364856781
//penguian_clean_dotplot,-0.8106666177952848
//eli_twisted_basepairs,0.21601775412216717
//aldo_loops_and_stacks,0.15067172153854896
//eli_direction_of_gc_pairs_in_multiloops_neckarea,0.09865127259968665
//eli_multiloop_similarity,0.00978877166914174
//eli_green_line,0.007680977135075769
//ding_quad_energy,0.00978877166906401
//berex_berex_loop_basic,0.019145505229272428
//eli_legal_placement_of_GUpairs,0.05120674277182907
//aldo_mismatch,-0.19581776453734906
//eli_red_line,-0.0018275413712541821
//eli_wrong_direction_of_gc_pairs_in_multiloops,-0.19731059143533014
//deivad_deivad_strategy,0.08548851359538262
//eli_direction_of_gc_pairs_in_multiloops,0.09865127290814124
//eli_no_blue_nucleotides_strategy,-0.8253549680053338
//berex_basic_test,-0.2751340813578167
//kkohli_test_by_kkohli,-0.06803891875360141
// => BadHelix3,0.10845716043121478
//SingleStrandTripleA,-0.2524611477466161
// => BadHairpin4Doublet,0.11621958146693934
// => BadHelix2,0.051953978499975795
//ViennaStructureComp,1.5272406678598347
//BadHairpin4Triplet,0.14470067211618373
// => BadSingleStrand3,0.1902229692459677
// => PctGCPairs,0.6973340034106996
// => Helix1,0.5804069633989485
// => JunctionClosing,-0.7514795405184205
//SingleStrand3,0.08819793678377201
//GRepeats,-0.05032464330299019
//ViennaMFENormalized,-0.49663024818195145
//GoodTetraloopSinglet,0.8372216739418228
//Helix4,0.6596505460982536
//BadJunction2_1_1,0.24122961401641185
//BadSingleStrand5,0.14447503330605593
//Helix3,0.29379550128941323
//GoodTetraloopDoublet,1.393824377152868
//GoodJunction2_1_1,-0.4598285526297501
//goodjunction3_1_1_1,-9.462788548128273
//GoodTetralooptriplet,1.916441326283273
//SingleStrand5,-0.11077156860161347
//NupackStructureComp,0.7965075352904527
//EliGreenLine,-0.004643042781698109
//GoodJunction3_3_3_3,1.1734967921515913
//Badjunction3_3_3_3,-0.097400093245612
//SpotStructureComp,0.4615538691139922
//BadSingleStrand4,0.05046284982294122
//Helix2,0.331594773993627
//GoodHairpin4,-1.0079078237755852
//AldoRepetition,0.15268020542080335
//AldoMismatch,0.2430050990096834
//RnastructureStructureComp,-0.23831987534823496
//NonCanonical,-0.0853402523094039
//BadHelix4,0.17028695543216954
//JunctionNonGCClosing,0.30753499424534314
