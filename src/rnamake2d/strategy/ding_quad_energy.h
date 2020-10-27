#ifndef RNAMAKE_DINGQUADENERGY_H
#define RNAMAKE_DINGQUADENERGY_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class DingQuadEnergy : public Strategy2D {
    private:

    public:
        DingQuadEnergy() : Strategy2D() {
            params_ = {-1.3405704811835677,
                       -0.7830146926406039,
                       -2.5760802053635405,
                       -11.287805083054593,
                       -1.3633074402117116,
                       -2.009150681671608,
                       7.852624100955059,
                       -3.1388969163494496,
                       -1.3541464922232402};
            name_ = "ding_quad_energy";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);

            auto GCEndingStems(0);
            auto AUEndingStems(0);
            auto TotalStems(0);

            for(auto ii = 0; ii < feature->elements.size(); ++ii) {
                const auto elem = feature->elements[ii];

                if(elem.type_ == RNAELEMENT::STACK) {
                    ++TotalStems;
                }

                for(auto ii = 0; ii < elem.quad_scores_.size(); ++ii) {
                    if(float(elem.quad_scores_[ii])/100.f > params_[0]) {
                        result += params_[1];
                    }
                    auto lengthOfStack = elem.get_stack_length();
                    for(auto ii = 0; ii < lengthOfStack - 1; ++ii) {
                        const auto pair = elem.get_pair_from_stack(ii, feature->sequence);
                        const auto nextPair = elem.get_pair_from_stack(ii + 1, feature->sequence);

                        if(pair == "GU" || pair == "UG" || nextPair == "GU" || nextPair == "UG") {
                            if(elem.quad_scores_.size() - 1 >= ii) {
                                if(float(elem.quad_scores_[ii]) / 100.f > params_[2] ) {
                                    result += params_[4];
                                }
                            }
                        }
                    }

                    if(lengthOfStack < params_[6]) {
                        for(auto ii = 0; ii < lengthOfStack -2 ; ++ii) {
                            const auto pair0 = elem.get_pair_from_stack(ii, feature->sequence);
                            const auto pair1 = elem.get_pair_from_stack(ii+1, feature->sequence);
                            const auto pair2 = elem.get_pair_from_stack(ii+2, feature->sequence);

                            if( (pair0 == "AU" || pair0 == "UA")
                            and (pair1 == "AU" || pair1 == "UA")
                            and (pair2 == "AU" || pair2 == "UA") ) {
                                result += params_[5];
                            }

                        }
                    }
                    const auto pairStart = elem.get_pair_from_stack(0, feature->sequence);
                    const auto pairEnd = elem.get_pair_from_stack(lengthOfStack - 1, feature->sequence);

                    if((pairStart == "GC" || pairStart == "CG" )
                    and (pairEnd == "GC" || pairEnd == "CG")) {
                        ++GCEndingStems;
                    }
                    if((pairStart == "AU" || pairStart == "UA")
                    and (pairEnd == "AU" || pairEnd == "UA")) {
                        ++AUEndingStems;
                    }
                }
            }
            const auto mutliplier = TotalStems - GCEndingStems;
            if(mutliplier > 0 && float(AUEndingStems) - params_[7] > 0.f) {
                result += params_[8]*(float(mutliplier) - std::min(float(AUEndingStems), params_[7]));
            }
            return result;
        }
    };

}
#endif // RNAMAKE_DINGQUADENERGY_H 