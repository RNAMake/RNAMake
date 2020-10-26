#ifndef RNAMAKE_GOODJUNCTION3_1_1_1_H
#define RNAMAKE_GOODJUNCTION3_1_1_1_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class GoodJunction3_1_1_1 : public Strategy2D {
    private:
        std::unordered_set<String> good_ = {
                            "GAG&CAG&CAC",
                            "GAG&CAC&GAC",
                            "CAG&CAC&GAG",
                            "GAC&GAC&GAC",
                            "CAC&GAC&GAG",
                            };
    public:
        GoodJunction3_1_1_1() : Strategy2D() {
            params_ = {9.18392777e+01, -4.13610363e-02};
            name_ = "GoodJunction3_1_1_1";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(params_[0]);

            for(const auto& motif : feature->motifs) {
               if(motif->mtype() != util::MotifType::NWAY) {
                   continue;
               }

               if(good_.find(motif->sequence) != good_.end()) {
                   result += params_[1];
               }
            }

            return result;
        }
    };

}
#endif // RNAMAKE_GOODJUNCTION3_1_1_1_H 