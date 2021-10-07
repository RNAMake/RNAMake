//headers for testing
#include "../common.hpp"
#include "../tools/motif_graph_builder.hpp"
#include "../tools/motif_tree_builder.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "base/backtrace.h"
#include "resources/resource_manager.h"
#include "secondary_structure/util.h"
#include "structure/is_equal.h"
#include "motif_data_structure/motif_graph.h"

#include <util/find_pair.h>

TEST_CASE( "Test Assembling Motifs together in Graph " ) {
    &HELIX.AVG.4&-4.6&0&9&ADE,A,1,A,,P 2.38596972398 -8.00026912158 -1.37124787419,OP2 3.46881960271 -7.05896897484 -1.72420519778,OP1 2.01703031066 -9.10959768342 -2.28908190787,O5' 1.05927229843 -7.18684286869 -1.0558599273,C5' -0.107434157857 -7.83538543977 -0.567513435256,C4' -1.21304329205 -6.83727493751 -0.341188918482,O4' -0.819597703626 -5.91016697849 0.701035000607,C3' -1.51076086087 -5.9452401513 -1.5351359488,O3' -2.37337794299 -6.56456651359 -2.46515296939,C1' -1.25143907704 -4.60526131267 0.379159396339,C2' -2.08871990659 -4.69509456902 -0.88851303494,O2' -3.43755164551 -4.91960919552 -0.507468699221,N1 0.567090643131 0.21096936527 0.18111835614,C2 -0.672007264061 -0.277626443567 0.304324114734,N3 -1.07266760087 -1.53987815685 0.34223741435,C4 -0.0256170205416 -2.3721305681 0.237468404972,C5 1.30695343645 -2.02045817283 0.102877858827,C6 1.60487173571 -0.642309208069 0.0780739489227,N6 2.83388527822 -0.130647558082 -0.0454416662578,N7 2.10099275925 -3.15915483768 0.0238431591396,C8 1.24024748773 -4.14772922951 0.101522911584,N9 -0.061566794199 -3.74786107642 0.22908689886,;URA,U,2,A,,P -1.61933461091 -6.32758894324 -4.40404052606,OP2 -0.240640796444 -6.16518264872 -4.93390994081,OP1 -2.58524059645 -7.18290848958 -5.14256087901,O5' -2.26019016605 -4.87891055562 -4.23398790423,C5' -3.67213563876 -4.70947472317 -4.03656630057,C4' -4.03159857262 -3.24128015665 -4.07620958631,O4' -3.32065453093 -2.54294723648 -3.01920405554,C3' -3.63776301755 -2.48682377928 -5.33611804295,O3' -4.56169508162 -2.67676004331 -6.40106248699,C1' -2.96636921007 -1.23947803859 -3.46053867428,C2' -3.59077295211 -1.04883698905 -4.8407658844,O2' -4.86421454221 -0.457742581499 -4.68825897078,N1 -1.49905015973 -1.14419423569 -3.50159738605,C2 -0.926774237848 0.118198987195 -3.51000928382,O2 -1.58460679946 1.14706440515 -3.5027502572,N3 0.447962930783 0.130343187871 -3.53592995566,C4 1.29380786314 -0.965762053566 -3.56427688999,O4 2.51198525161 -0.793287480544 -3.57481846959,C5 0.62730214243 -2.23209970501 -3.56573634194,C6 -0.709547577961 -2.27867135818 -3.53248541544,;URA,U,3,A,,P -4.2165942583 -2.32319938323 -8.29394483185,OP2 -2.94974970513 -3.06154754915 -8.48986301551,OP1 -5.35774789579 -2.54287029393 -9.21241255238,O5' -3.85997977877 -0.764351010473 -8.2523473754,C5' -4.897254748 0.203400688532 -8.28903414556,C4' -4.30108831532 1.58700840747 -8.18664473963,O4' -3.51419681588 1.67979126766 -6.97023298965,C3' -3.30207834065 1.922813164 -9.28369807895,O3' -3.9641878966 2.31361211714 -10.4772242156,C1' -2.40446725626 2.53192084636 -7.19840065041,C2' -2.53337572327 3.06007020859 -8.62683822969,O2' -3.24060519448 4.28129893227 -8.63157835554,N1 -1.12371841986 1.7817695265 -6.97215873391,C2 0.0264273703221 2.50381589233 -6.71776686964,O2 0.0493893211892 3.72072581343 -6.66194681547,N3 1.15612592084 1.74685244475 -6.52763712066,C4 1.24882858032 0.36678425976 -6.56916686177,O4 2.32918446614 -0.17266002565 -6.3729760633,C5 0.00852200312181 -0.321529369711 -6.84590971749,C6 -1.10423547412 0.399363312414 -7.03147951798,;GUA,G,4,A,,P -3.31822283877 2.00140392187 -11.9114176766,OP2 -2.78183052197 0.624550668562 -11.9321211852,OP1 -4.31221514664 2.41152616892 -12.9303283013,O5' -2.07300822276 3.00699424535 -11.9616922146,C5' -2.27997758731 4.40567265631 -12.1769839231,C4' -0.981753769229 5.18555245058 -12.0550595249,O4' -0.455333191491 5.0900900384 -10.7004529242,C3' 0.191044909974 4.6966179921 -12.9026470294,O3' 0.0484368902736 5.0125668221 -14.2919464975,C1' 0.959850702052 5.20848312898 -10.7567321805,C2' 1.31308600457 5.47161050737 -12.2230251348,O2' 1.28299766212 6.85965223516 -12.499214454,N1 5.32395030992 3.14930017638 -9.16003760091,C2 5.05833494743 4.44279897514 -9.54113097111,N2 6.08398053361 5.30212500277 -9.49348975051,N3 3.86477379873 4.86101304896 -9.93731432283,C4 2.94003960702 3.86317816608 -9.92623523959,C5 3.10949569255 2.54318800002 -9.56653342644,C6 4.39385665808 2.11120514817 -9.14150891465,O6 4.75399966679 0.982471023003 -8.78058290433,N7 1.91166867176 1.84388986134 -9.69463577576,C8 1.05871033337 2.73768962367 -10.1207505591,N9 1.60862251341 3.98761966582 -10.2763328694,;:CYT,C,5,A,,P 14.0752573663 1.19626348221 -4.57140079991,OP2 12.9949638704 0.330912002986 -4.04672283309,OP1 15.1354368515 1.6914232433 -3.67116200874,O5' 13.3312482803 2.41089749464 -5.31098149311,C5' 14.0570401658 3.47182417062 -5.94420285479,C4' 13.1088342879 4.43938449053 -6.63055838059,O4' 12.4334604786 3.74022226474 -7.70323342041,C3' 11.9505142274 4.95192917865 -5.7827204847,O3' 12.3379854394 6.03669884067 -4.96537638822,C1' 11.1377543306 4.2858630117 -7.88776718182,C2' 10.9683222226 5.39238206517 -6.85798070182,O2' 11.2984103272 6.64082196413 -7.4353760992,N1 10.0978697976 3.21290540438 -7.76638031024,C2 8.82122402458 3.43051938761 -8.31743642865,O2 8.5752749841 4.49601002993 -8.88970134739,N3 7.87999802282 2.45957345199 -8.20771736105,C4 8.17242744926 1.32139310906 -7.57991237347,N4 7.21294119429 0.403856225057 -7.51317313414,C5 9.46093342163 1.07332915105 -7.00828633268,C6 10.3826796467 2.03651483174 -7.12509573486,;ADE,A,6,A,,P 11.5029614721 6.37921491207 -3.64380532962,OP2 11.1891715519 5.11876412592 -2.92287980029,OP1 12.228934435 7.4622809749 -2.94966651153,O5' 10.1211607291 6.93299744635 -4.22403023034,C5' 10.0193471574 8.24921009723 -4.75911112748,C4' 8.58586819007 8.53194034072 -5.15628891899,O4' 8.17572755955 7.59695201847 -6.18649213465,C3' 7.55480640286 8.31254509991 -4.06279509699,O3' 7.52150636112 9.40927742302 -3.16951679094,C1' 6.78334980448 7.35209878696 -6.0669735708,C2' 6.29132991422 8.20320161133 -4.90362364626,O2' 5.84495657434 9.4529482745 -5.379269895,N1 3.43547638441 3.42536737617 -6.33578181017,C2 3.30832953586 4.72860762965 -6.62178316542,N3 4.22180277767 5.69702152019 -6.53131264148,C4 5.39427414687 5.21557665669 -6.0785800905,C5 5.67876043561 3.90179850389 -5.74050661853,C6 4.62786744226 2.96513284734 -5.88504881089,N6 4.74839421627 1.66362246115 -5.59569753442,N7 6.98860635588 3.78102602179 -5.31072926115,C8 7.46495815174 4.99728435949 -5.40034286979,N9 6.56140987772 5.9163535804 -5.85995815868,;ADE,A,7,A,,P 6.87130720721 9.08500310426 -1.42946847985,OP2 7.33957109523 7.77076651482 -0.919253440677,OP1 7.13910652139 10.3084410351 -0.635046948371,O5' 5.29689203692 8.98682977991 -1.66247128938,C5' 4.54552090171 10.1167117278 -2.15864276339,C4' 3.13514257039 9.69454197424 -2.5099958696,O4' 3.17201318475 8.73646130384 -3.60277301677,C3' 2.36704900488 8.9655009794 -1.42138244516,O3' 1.8160914052 9.83791436539 -0.44614617467,C1' 2.12958608927 7.78825333186 -3.44087518299,C2' 1.30914209 8.22620049109 -2.22615315718,O2' 0.231686505363 9.05115646313 -2.62742163156,N1 1.46383980697 2.69250393547 -3.5978152523,C2 0.645989307102 3.7180530853 -3.87711849072,N3 0.882232670639 5.02129641499 -3.82085476248,C4 2.1402061796 5.2553032543 -3.40998529318,C5 3.09658668606 4.31007591569 -3.08003947446,C6 2.72222982271 2.96040080944 -3.19044550612,N6 3.53838242933 1.93721171506 -2.92147227259,N7 4.28215151219 4.92253180913 -2.69573799121,C8 4.02382255306 6.20277854697 -2.80457136543,N9 2.74433411858 6.47698148068 -3.23008707099,;URA,U,8,A,,P 2.23462174995 9.75765908961 1.05069203841,OP2 3.44418940313 8.97369564525 1.40018268874,OP1 1.99997118521 11.0293394927 1.78127892701,O5' 0.935166983354 8.82861736902 1.16715642995,C5' -0.352932549123 9.30719969143 0.800852238578,C4' -1.35799739772 8.18926490797 0.638812237486,O4' -0.953659205185 7.31117003505 -0.443793129497,C3' -1.52969869889 7.24063274032 1.81669924849,O3' -2.32331164609 7.75545163342 2.86830226444,C1' -1.35556092296 5.98026131267 -0.157159396339,C2' -2.12899976502 6.00949843406 1.15631101813,O2' -3.49972891075 6.2168188897 0.858125501297,N1 -0.153213789785 5.1236696498 -0.0448602279318,C2 -0.3361466684 3.74984680465 -0.0550713749148,O2 -1.4327067546 3.22355011587 -0.161890523941,N3 0.82064944115 3.01554437464 0.0620327343649,C4 2.10248374306 3.51499390523 0.184137260446,O4 3.0558347086 2.74442837819 0.281519433718,C5 2.20336076403 4.94521546317 0.191157110641,C6 1.09462061766 5.67621861075 0.0804589472231,;:&A2-A7,10.6310593834 -1.59552712743 -13.879367919;0.855774221129 -0.501127336115 0.128397963278 0.496017150236 0.865297359955 0.0719668594572 -0.147192870846 0.00213973445563 0.989103958179 ;10.4302174859 3.01608691253 -16.4422121526 13.6849137734 -5.43479406756 -11.3965562488 ,cW-W@A1-A8,-18.5831468979 -0.387081470873 0.202736602578;1.00002299863 9.99985999616e-05 -1.000031985e-08 -8.39964971388e-05 0.999985999616 -0.000100003199117 1.30088343956e-05 0.000107002088103 1.00003200117 ;-1.25143907704 -4.60526131267 0.379159396339 -1.35556092296 5.98026131267 -0.157159396339 ,cW-W@A3-A6,10.8566629708 -5.03063686449 -15.6265717355;0.424928846447 -0.88913356166 0.169998832932 0.872220159896 0.452432466292 0.18601982 -0.24226164636 0.0692813842788 0.967755165276 ;13.5678539381 -1.03346785891 -17.9479143658 10.854642924 -9.74160850657 -12.8747499281 ,cW-W@A5-A4,9.21317219454 -8.28261135837 -17.2980176363;-0.141821813685 -0.948090545269 0.284603433197 0.951852600552 -0.0516761883795 0.302121902079 -0.271682747805 0.31377837676 0.909848709208 ;6.13847433011 -12.0877438229 -14.857508827 14.0393345317 -6.45908704795 -19.1663997865 ,cW-W@&1 3 &AUUG_LLLL_CAAU_RRRR CAAU_LLLL_AUUG_RRRR &99!assembled!assembled!A,(,1,A,;U,(,2,A,;U,(,3,A,;G,(,4,A,;|C,),5,A,;A,),6,A,;A,),7,A,;U,),8,A,;|!4 3@2 5@0 7@1 6@!2 0 !AUUG_LLLL_CAAU_RRRR CAAU_LLLL_AUUG_RRRR &&

    SUBCASE("test adding motifs") {
        auto mg = motif_data_structure::MotifGraph();
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL");
        auto m3 = resources::Manager::instance().motif("HELIX.IDEAL");
        
        SUBCASE("cannot find end if there is not parent") {
            REQUIRE_THROWS_AS(mg.add_motif(m1, -1, "A1-A8"), motif_data_structure::MotifGraphException);
        }
        
        mg.add_motif(m1);

        CHECK(mg.size() == 1);
        
        // can never use parent_end_index=0 for a graph as that is where that node
        // is already connected to another node
        REQUIRE_THROWS_AS(mg.add_motif(m2, -1, 0), motif_data_structure::MotifGraphException);
        
        // catches invalid parent_index
        REQUIRE_THROWS_AS(mg.add_motif(m2, 2), motif_data_structure::MotifGraphException);
        
        // invalid parent_end_index, has only 0 and 1
        REQUIRE_THROWS_AS(mg.add_motif(m2, -1, 3), motif_data_structure::MotifGraphException);
        
        // invalid parent_end_name, is the name of end 0
        REQUIRE_THROWS_AS(mg.add_motif(m2, -1, "A4-A5"), motif_data_structure::MotifGraphException);
        
        // invalid parent_end_name, cannot be found as an end in motif
        REQUIRE_THROWS_AS(mg.add_motif(m2, -1, "FAKE"), motif_data_structure::MotifGraphException);
    }
    
    SUBCASE("test removing motifs") {
        auto mg = motif_data_structure::MotifGraph();
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m3 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        
        mg.add_motif(m1);
        mg.add_motif(m2);
        
        mg.remove_motif(1);
        CHECK(mg.size() == 1);
        
        mg.add_motif(m2);
        mg.add_motif(m3);
        mg.remove_level(0);
        CHECK(mg.size() == 0);
        
    }
    
    SUBCASE("test getting nodes") {
        auto mg = motif_data_structure::MotifGraph();
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        mg.add_motif(m1);
        mg.add_motif(m2);
        
        REQUIRE_NOTHROW(mg.get_node(0));
        REQUIRE_THROWS(mg.get_node(10));
        
    }
    
    SUBCASE("test stringifying the motif graph yeilds indentical motifs when reloaded") {
        auto builder = MotifGraphBuilder();
        auto mg = builder.build(5);
        auto s = mg->to_str();
        
        auto mg2 = std::make_shared<motif_data_structure::MotifGraph>(s, motif_data_structure::MotifGraphStringType::MG);
        CHECK(mg->size() == mg2->size());

        for(int i = 0; i < mg->size(); i++) {
            auto atoms1 = mg->get_node(i)->data()->atoms();
            auto atoms2 = mg2->get_node(i)->data()->atoms();
            CHECK(are_atom_vectors_equal(atoms1, atoms2));
            
        }
        
        auto struc = mg2->get_structure();
        CHECK(struc->chains().size() == 2);
    }
    
    SUBCASE("test compatibility with python stringification and topology serialization") {
        auto path = base::base_dir() + "/unittests/unittest_resources/motif_graph/";
        auto lines =base::get_lines_from_file(path + "test.mg");
        auto mg = motif_data_structure::MotifGraph(lines[0], motif_data_structure::MotifGraphStringType::MG);
        auto s = mg.get_structure();
        
        lines =base::get_lines_from_file(path + "base_mg.mg");
        mg = motif_data_structure::MotifGraph(lines[0], motif_data_structure::MotifGraphStringType::MG);
        s = mg.get_structure();
        
    }

    SUBCASE("test copying motif graph") {
        auto builder = MotifGraphBuilder();
        auto mg = builder.build(5);
        
        auto mg_copy = std::make_shared<motif_data_structure::MotifGraph>(*mg);
        CHECK(mg->size() == mg_copy->size());
        
        for(int i = 0; i < mg->size(); i++) {
            auto atoms1 = mg->get_node(i)->data()->atoms();
            auto atoms2 = mg_copy->get_node(i)->data()->atoms();
            CHECK(are_atom_vectors_equal(atoms1, atoms2));
            
        }
        for(auto const & n : *mg_copy) {
            int count = 0;
            for(auto const & c : n->connections()) {
                if(c != nullptr) { count++; }
            }
            CHECK(count != 0);
        }
        
        auto rna_struct = mg_copy->get_structure();
        CHECK(rna_struct->chains().size() == 2);
        
        mg->replace_ideal_helices();
        auto mg_copy_2 = std::make_shared<motif_data_structure::MotifGraph>(*mg);
        auto dss = mg_copy_2->designable_secondary_structure();
        secondary_structure::fill_basepairs_in_ss(dss);
        mg_copy_2->replace_helical_sequence(dss);

        
        
    }
 
    SUBCASE("test replacing idealized helices") {
        auto mg = motif_data_structure::MotifGraph();
        auto m = resources::Manager::instance().motif("HELIX.IDEAL.6");
        mg.add_motif(m);
        mg.replace_ideal_helices();
        
        CHECK(mg.size() == 7);
    
        auto builder = MotifGraphBuilder();
        auto mg2 = builder.build(2);
        
        SUBCASE("make sure complex build is producing the right number of motifs") {
        
            auto expected_total = 0;
            for(auto const & n : *mg2) {
                if(n->data()->mtype() != util::MotifType::HELIX) {
                    expected_total += 1;
                }
                else {
                    expected_total += n->data()->basepairs().size() -1;
                }
            }
        
            mg2->replace_ideal_helices();
            CHECK(mg2->size() == expected_total);
        }
        
        auto s = mg2->get_structure();
        
        CHECK(s->chains().size() == 2);
    }
  
    SUBCASE("test replacing idealized helices 2") {
        auto mg = std::make_shared<motif_data_structure::MotifGraph>();
        auto m = resources::Manager::instance().motif("HELIX.IDEAL.6");
        m->move(math::Point(40, 0, 0));
        mg->add_motif(m);
        
        auto new_mg = std::make_shared<motif_data_structure::MotifGraph>(*mg);
        new_mg->replace_ideal_helices();
        
        auto struc1 = mg->get_structure();
        auto struc2 = new_mg->get_structure();
        
        auto d1 = struc1->ends()[1]->d();
        auto ds_2 = math::Points{struc1->ends()[0]->d(), struc1->ends()[1]->d()};
        
        auto dist_1 = d1.distance(ds_2[0]);
        auto dist_2 = d1.distance(ds_2[1]);
        auto min = dist_1 < dist_2 ? dist_1 : dist_2;
        
        CHECK(min < 0.1);
    }
    
    SUBCASE("test replacing helices with new sequence") {
        auto mg = motif_data_structure::MotifGraph();
        auto m = resources::Manager::instance().motif("HELIX.IDEAL.6");
        mg.add_motif(m);
        mg.replace_ideal_helices();

        auto dss = mg.designable_secondary_structure();
        secondary_structure::fill_basepairs_in_ss(dss);
        
        REQUIRE_NOTHROW(mg.replace_helical_sequence(dss));
        
        auto builder = MotifGraphBuilder();
        auto mg2 = builder.build(2);
        mg2->replace_ideal_helices();

        dss = mg2->designable_secondary_structure();
        
        secondary_structure::fill_basepairs_in_ss(dss);

        REQUIRE_NOTHROW(mg2->replace_helical_sequence(dss));
        CHECK(mg2->sequence() == dss->sequence());
        
    }
    
    SUBCASE("test replacing helices with new sequence 2") {
        auto mg = std::make_shared<motif_data_structure::MotifGraph>();
        auto m = resources::Manager::instance().motif("HELIX.IDEAL.6");
        m->move(math::Point(40, 0, 0));
        mg->add_motif(m);
        mg->replace_ideal_helices();
        
        auto new_mg = std::make_shared<motif_data_structure::MotifGraph>(*mg);
        auto dss = new_mg->designable_secondary_structure();
        secondary_structure::fill_basepairs_in_ss(dss);
        new_mg->replace_helical_sequence(dss);

        auto struc1 = mg->get_structure();
        auto struc2 = new_mg->get_structure();
        
        auto d1 = struc1->ends()[1]->d();
        auto ds_2 = math::Points{struc1->ends()[0]->d(), struc1->ends()[1]->d()};
        
        auto dist_1 = d1.distance(ds_2[0]);
        auto dist_2 = d1.distance(ds_2[1]);
        auto min = dist_1 < dist_2 ? dist_1 : dist_2;
        
        CHECK(min < 10);
    }
    
    SUBCASE("test get end for easy building ") {
        auto base_path = base::base_dir() + "/apps/simulate_tectos/resources/";
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL.2");

        resources::Manager::instance().add_motif(base_path+"GAAA_tetraloop", "ttr");
        
        auto mg = motif_data_structure::MotifGraph();
        auto ttr_m = resources::Manager::instance().motif("ttr", "", "A229-A245");
        mg.add_motif(m1);
        mg.add_motif(ttr_m);
        mg.add_motif(m2);
        
        CHECK(mg.size() == 3);
        
        REQUIRE_NOTHROW(mg.get_available_end(1, "A222-A251"));

        REQUIRE_THROWS_AS(mg.get_available_end(0), motif_data_structure::MotifGraphException);
        REQUIRE_THROWS_AS(mg.get_available_end(2, "A4-A5"), motif_data_structure::MotifGraphException);

        SUBCASE("try getting end by the name of the motif and end name") {
        
            REQUIRE_NOTHROW(mg.get_available_end("ttr", "A222-A251"));
            REQUIRE_THROWS_AS(mg.get_available_end("ttr", "FAKE_END"), motif_data_structure::MotifGraphException);
            REQUIRE_THROWS_AS(mg.get_available_end("FAKE_MOTIF", "A222-A251"), motif_data_structure::MotifGraphException);
            REQUIRE_THROWS_AS(mg.get_available_end("HELIX.IDEAL.2", "A1-A8"), motif_data_structure::MotifGraphException);
            REQUIRE_THROWS_AS(mg.get_available_end("ttr", "A149-A154"), motif_data_structure::MotifGraphException);

        }
        
    }
    
    SUBCASE("test connecting motifs") {
        auto mg = motif_data_structure::MotifGraph();
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m3 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto nway = resources::Manager::instance().motif("NWAY.1GID.0");
        mg.add_motif(m1);
        mg.add_motif(nway);
        mg.add_motif(m2);

        // try connecting through 0th end position
        REQUIRE_THROWS_AS(mg.add_connection(1, 2, "A138-A180", ""), motif_data_structure::MotifGraphException);
        
        // try connecting thru an already used end position
        REQUIRE_THROWS_AS(mg.add_connection(1, 2, "A141-A162", ""), motif_data_structure::MotifGraphException);
        
        mg.add_connection(1, 2, "", "");
        auto s = mg.get_structure();
        CHECK(s->chains().size() == 1);

        REQUIRE_THROWS_AS(mg.add_motif(m3, -1, 1), motif_data_structure::MotifGraphException);
        CHECK(mg.add_motif(m3) == -1);
        
        REQUIRE_THROWS_AS(mg.add_connection(1, 2, "", ""), motif_data_structure::MotifGraphException);
        REQUIRE_THROWS_AS(mg.add_connection(1, 2, "", m1->ends()[0]->name()), motif_data_structure::MotifGraphException);
        
    }
    
    SUBCASE("test adding motif tree to motif graph") {
        
        auto g = std::make_shared<BuilderGraph>();
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        g->add_node(util::MotifType::HELIX, 2);
        g->add_node(util::MotifType::NWAY, 3);
        g->add_node(util::MotifType::HELIX, 2, 1, 1);
        g->add_node(util::MotifType::HELIX, 2, 1, 2);
        
        auto builder = MotifTreeBuilder(g);
        auto mt = builder.build();
        
        auto mg = motif_data_structure::MotifGraph();
        mg.add_motif(m1);
        REQUIRE_NOTHROW(mg.add_motif_tree(mt));
        CHECK(mg.size() == 5);
        CHECK(mg.get_node(2)->data()->mtype() == util::MotifType::NWAY);
        
        auto mg2 = motif_data_structure::MotifGraph();
        mg2.add_motif(m1);
        REQUIRE_NOTHROW(mg2.add_motif_tree(mt, 0, "A1-A8"));
        CHECK(mg2.size() == 5);

        
    }
}



















