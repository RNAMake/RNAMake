#ifndef RNAMAKE_DESIGN_H
#define RNAMAKE_DESIGN_H

#include <iostream>
#include <utility>
#include <vector>
#include <sstream>

#include <base/types.h>
#include <rnamake2d/SSMotif.h>
#include <rnamake2d/feature_generator2d.h>
#include <eternabot/strategy/a_basic_test.h>
#include <rnamake2d/rna_element.h>
#include <rnamake2d/rnamake2d.fwd.hh>

namespace rnamake2d {


    class Design {
        // what do we want to be doing here??
        FeatureGenerator2d generator;
        Feature2DOP best;
    public:
        Feature2DOP mutant_;
    private:
        secondary_structure::PoseOP best_pose;
        secondary_structure::PoseOP mutant_pose;
        secondary_structure::Parser parser;
        Reals scores;
        static int design_ct;
        int num_moves{0};
        float bp_score_{0.f};
        double score_{0.};
        const String target_;
        String sequence_; // the current sequence
        String candiate_; // the candidate
        String constraint_;
        const Ints pairmap;
        const int id;
        friend NemoSampler;
    public:
        explicit
        Design(String targ, String seq = "") : id(design_ct++),
                                                 target_(targ),
                                                 sequence_(seq),
                                                 best(nullptr),
                                                 mutant_(nullptr),
                                                 generator(),
                                                 pairmap(get_pairmap_from_secstruct(targ))
                                                 {
            while(target_.size() - sequence_.size()) {
                 sequence_.push_back('N');
            }
            constraint_ = sequence_;
        }

    public:
        void
        initialize_features() ;

    public:
        void
        initialize_mutant(String& ) ;

    public:
        void
        accept(double mutant_score);

    public:
        void
        accept_bp(double mutant_score) ;

    public:
        float
        score() const {
            return score_;
        }
    public:
        auto
        bp_score() const {
            return bp_score_;
        }

    public:
        void
        bp_score(float score) {
            bp_score_ = score;
        }

    public:
        void
        reject() {
            scores.push_back(score_);
            ++num_moves;
            candiate_ = sequence_;
        }

    public:
        const Feature2DOP
        mutant() const {
            return mutant_;
        }

    public:
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

    public:
        String
        target() const {
            return target_;
        }
        String
        candidate() const {
            return candiate_;
        }

        String
        sequence() const {
            return sequence_;
        }
        String
        to_str() const {
            auto ss = std::stringstream{};
            ss<<id<<","<<target_<<","<<sequence_<<","<<score_<<","<<num_moves<<"\n";
            return ss.str();
        }
   };


}

using Designs = std::vector<rnamake2d::Design>;

namespace rnamake2d {
    Designs
    make_designs(String target, String  sequence="",int number=1);
}

#endif  // RNAMAKE_DESIGN_H