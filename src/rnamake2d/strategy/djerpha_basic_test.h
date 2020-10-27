#ifndef RNAMAKE_DJERPHABASICTEST_H
#define RNAMAKE_DJERPHABASICTEST_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class DjerphaBasicTest : public Strategy2D {
    private:

    public:
        DjerphaBasicTest() : Strategy2D() {
            params_ = {0.5167815695227542,
                       29.124574301094423,
                       -1.4576430236665883,
                       0.008240741000965893,
                       5.356389162410184,
                       47.140830683223825,
                       200.69890875030296 };
            name_ = "djerpha_basic_test";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto score(100.f);
            const auto total_pairs = feature->gu + feature->gc + feature->ua;

            if(total_pairs > 0.f) {
                score -= std::abs(feature->ua / total_pairs - params_[0])*params_[1];
            }
            const auto target_fe = params_[2]*total_pairs;

            score -= std::abs(target_fe - feature->fe) *params_[3];

            if(feature->meltpoint < params_[5]) {
                score -= std::abs(feature->meltpoint - params_[5]) * params_[4];
            } else if (feature->meltpoint > params_[6]) {
                score -= std::abs(feature->meltpoint - params_[6]) * params_[4];
            }

            return score;
        }
    };

}
#endif // RNAMAKE_DJERPHABASICTEST_H 