#ifndef RNAMAKE_ELINOBLUENUCLEOTIDESSTRATEGY_H
#define RNAMAKE_ELINOBLUENUCLEOTIDESSTRATEGY_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class EliNoBlueNucleotidesStrategy : public Strategy2D {
    private:

    public:
        EliNoBlueNucleotidesStrategy() : Strategy2D() {
            params_ = {0.40219824218749967};
            name_ = "eli_no_blue_nucleotides_strategy";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);

            for(auto ii = 0; ii < feature->elements.size(); ++ii) {
                const auto &elem = feature->elements[ii];

                if (elem.type_ == RNAELEMENT::LOOP) {
                    const auto loop_groups = elem.get_loop_groups();
                    if (loop_groups.size() > 2 ||
                        (loop_groups.size() == 2 && elem.parent_ == nullptr)) {
                        for (auto kk = 0; kk < loop_groups.size(); ++kk) {
                            for (auto mm = 0; mm < loop_groups[kk].size(); ++mm) {
                                if (feature->sequence[loop_groups[kk][mm]] == 'U') {
                                    result -= params_[0];
                                }
                            }
                        }
                    }
                }
            }
            return result;
        }
    };

}
#endif // RNAMAKE_ELINOBLUENUCLEOTIDESSTRATEGY_H 