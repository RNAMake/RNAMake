//
// Created by cjurich on 10/19/20.
//

#ifndef RNAMAKE_BAD_HELIX2_H
#define RNAMAKE_BAD_HELIX2_H

#include <unordered_set>

#include <rnamake2d/strategy2d.h>

namespace rnamake2d {

    class BadHelix2 : public Strategy2D {
        std::unordered_set<String> bad_ = {
                "AU&AU",
                "CA&UG",
                "GC&GU",
                "GU&GC",
                "UU&AA"
        };
        float penalty_ = 9.88953488;
    public:
        BadHelix2() : Strategy2D() {

            name_ = "BadHelix2";
        }

    public:
        float
        score(Feature2DOP const& feature) override {
            auto helix_ct(0);
            for(auto& motif : feature->motifs) {
                if(motif->token() != "Helix2") {
                    continue;
                }
                if(bad_.find(motif->sequence) != bad_.end()) {
                    ++helix_ct;
                }
            }

            return 100.f - penalty_*helix_ct;
        }


    };


}

#endif //RNAMAKE_BAD_HELIX2_H
