#ifndef RNAMAKE_BAD_JUNCTION2_1_1_H
#define RNAMAKE_BAD_JUNCTION2_1_1_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class BadJunction2_1_1 : public Strategy2D {
    private:
        std::unordered_set<String> bad_ = {
                "CGC&GGG",
                "CGU&AGG",
                "GGU&AGC",
                "GGG&CGC",
                "UGG&CGA",
                "AGG&CGU",
                "CGA&UGG",
                "CAC&GAG",
                "CGG&CGG",
                "GGC&GGC",
                "AGA&UGU",
                "CUG&CUG",
                "UGA&UGA",
                "AGU&AGU",
                "GCU&AUC",
                "UGU&AGA",
                "UUU&AUG",
                "GUU&AUU",
                "UAA&UAA",
                "CAU&AAG",
                "ACA&UCU",
        };
    public:
        BadJunction2_1_1() : Strategy2D() {
            name_ = "BadJunction2_1_1";
            params_ =  {4.81526177};
        }

    public:
        float
        score( Feature2DOP const & feature) override {
            auto result(100.f);

            for(const auto& motif : feature->motifs) {
                if(motif->mtype() != util::MotifType::NWAY) {
                    continue;
                }

                if(bad_.find(motif->sequence) != bad_.end()) {
                    result -= params_[0];
                }
            }

            return result;
        }
    };
}

#endif // RNAMAKE_BAD_JUNCTION2_1_1_H