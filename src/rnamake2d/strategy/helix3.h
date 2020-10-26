#ifndef RNAMAKE_HELIX_3_H
#define RNAMAKE_HELIX_3_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class Helix3 : public Strategy2D  {
    private:
        std::unordered_set<String> good_ = {
                "CCA&UGG",
                "CCC&GGG",
                "CUG&CAG",
                "CUG&CGG",
                "GUG&CAC",
        };

    public:
        Helix3() : Strategy2D() {
            name_ = "Helix3";
            params_ = {92.27963326,  0.23109286,  0.80329228};

        }

    public:
        float
        score( Feature2DOP const & feature ) override {
            auto result(params_[0]);

            for(const auto& motif : feature->motifs) {
                if(motif->mtype() != util::MotifType::HAIRPIN || motif->buffer_ != 3) {
                    continue;
                }

                if(good_.find(motif->sequence) != good_.end()) {
                    result += params_[1] ;
                } else {
                    result -= params_[2];
                }
            }

            return result;
        }
    };
}

#endif // RNAMAKE_HELIX_3_H