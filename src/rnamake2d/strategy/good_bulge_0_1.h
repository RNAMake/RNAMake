#ifndef RNAMAKE_GOODBULGE_0_1_H
#define RNAMAKE_GOODBULGE_0_1_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class GoodBulge_0_1 : public Strategy2D {
    private:
        std::unordered_set<String> good_ = {
                "GG&CAC",
                "CG&CAG",
        };
    public:
        GoodBulge_0_1() : Strategy2D() {
            params_ = {80., 1.};
            name_ = "GoodBulge_0_1";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(params_[0]);

            for(const auto& junc : feature->motifs) {
                if(junc->mtype() != util::MotifType::NWAY &&
                    junc->mtype() != util::MotifType::TWOWAY
                ) {
                    continue;
                }
                if(good_.find(junc->sequence) != good_.end() ) {
                    result += params_[1];
                }
            }

            return result;
        }
    };

}
#endif // RNAMAKE_GOODBULGE_0_1_H 