#ifndef __RNAMAKE_DESIGN_H__
#define __RNAMAKE_DESIGN_H__

#include <iostream>
#include <utility>
#include <vector>
#include <sstream>

#include <base/types.h>
#include <rnamake2d/SSMotif.h>
#include <rnamake2d/feature_generator2d.h>
#include <eternabot/strategy/a_basic_test.h>

namespace rnamake2d {
    struct Design {
        FeatureGenerator2d generator;
        Feature2DOP best;
        Feature2DOP mutant;
        secondary_structure::PoseOP best_pose;
        secondary_structure::PoseOP mutant_pose;
        secondary_structure::Parser parser;
        Reals scores;
        static int design_ct;
        int num_moves{0};

        double score{0.};
        const String target;
        String sequence; // the current sequence
        String candiate; // the candidate
        const int id;

        explicit
        Design(String targ, String seq = "") : id(design_ct++),
                                                 target(targ),
                                                 sequence(seq),
                                                 best(nullptr),
                                                 mutant(nullptr),
                                                 generator()
                                                 {
            while(target.size() - sequence.size()) {
                 sequence.push_back('N');
            }
        }

        void
        initialize_features() ;

        void
        initialize_mutant() ;

        void
        accept(double mutant_score) {
            score = mutant_score;
            scores.push_back(mutant_score);
            ++num_moves;
            sequence = candiate;
            best = mutant;
            mutant = nullptr;
            best_pose = mutant_pose;
            mutant_pose = nullptr;
        }

        void
        reject() {

            scores.push_back(score);
            ++num_moves;
            candiate = sequence;

        }

        void
        update(bool mutate=false) {
//            for(auto& m : motifs) {
//                    m->full_mutated_sequence_ = candiate;
//                    m->full_sequence_ = sequence;
//            }
//
//            for(auto& m : motifs) {
//                m->build_sequence(mutate);
//            }
//            for(const auto& m : motifs) {
//                m->show();
//            }
        }

        String
        to_str() const {
            auto ss = std::stringstream{};
            ss<<id<<","<<target<<","<<sequence<<","<<score<<","<<num_moves<<"\n";
            return ss.str();
        }
   };

    struct Stats {
        int gc;
        int gu;
        int au;
    };

}

using Designs = std::vector<rnamake2d::Design>;

namespace rnamake2d {
    Designs
    make_designs(String target, String  sequence="",int number=1);
}

#endif  // __RNAMAKE_DESIGN_H__