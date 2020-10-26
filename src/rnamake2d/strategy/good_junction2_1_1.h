#ifndef RNAMAKE_GOOD_JUNCTION2_1_1_H
#define RNAMAKE_GOOD_JUNCTION2_1_1_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class GoodJunction2_1_1 : public Strategy2D {
    private:
        std::unordered_set<String> good_ = {
                                "CAG&CAG",
                                "GGA&UGC",
                                "GAC&GAC",
                            };

    public:
        GoodJunction2_1_1() : Strategy2D() {
            params_ = {91.87936627, -2.05195961};
            name_ = "GoodJunction2_1_1";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(params_[0]);
            for(const auto& motif : feature->motifs) {
                if(motif->mtype() == util::MotifType::NWAY &&
                    good_.find(motif->sequence) != good_.end()
                )  {
                    result += params_[1];
                }
            }
            return result;
        }
    };

}

#endif // RNAMAKE_GOOD_JUNCTION2_1_1_H