#ifndef RNAMAKE_QUASISPECIESTESTBYREGIONBOUNDARIES_H
#define RNAMAKE_QUASISPECIESTESTBYREGIONBOUNDARIES_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class QuasispeciesTestByRegionBoundaries : public Strategy2D {
    private:

    public:
        QuasispeciesTestByRegionBoundaries() : Strategy2D() {
            params_ = {1};
            name_ = "quasispecies_test_by_region_boundaries";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);

            for (auto ii = 0; ii < feature->elements.size(); ++ii) {
                if(feature->elements[ii].type_ == RNAELEMENT::STACK) {
                    if(feature->elements[ii].indices_.size() >= 4) {
                        if(feature->elements[ii].get_pair_from_stack(0, feature->sequence) ==
                        feature->elements[ii].get_pair_from_stack(1, feature->sequence)
                        ) {
                            result -= params_[0];
                        }
                    }
                    if(feature->elements[ii].indices_.size() >= 6) {
                        const auto last_pair = feature->elements[ii].indices_.size() / 2 - 1;
                        if(feature->elements[ii].get_pair_from_stack(last_pair, feature->sequence)
                        == feature->elements[ii].get_pair_from_stack(last_pair - 1, feature->sequence)) {
                            result -= params_[0];
                        }
                    }
                }
            }

            return result;
        }
    };

}
#endif // RNAMAKE_QUASISPECIESTESTBYREGIONBOUNDARIES_H 