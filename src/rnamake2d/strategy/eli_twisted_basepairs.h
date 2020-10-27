#ifndef RNAMAKE_ELITWISTEDBASEPAIRS_H
#define RNAMAKE_ELITWISTEDBASEPAIRS_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class EliTwistedBasepairs : public Strategy2D {
    private:

    public:
        EliTwistedBasepairs() : Strategy2D() {
            params_ = {};
            name_ = "eli_twisted_basepairs";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto score(100.f);

            for(auto ii = 0 ; ii < feature->elements.size(); ++ii) {
                if(feature->elements[ii].type_ == RNAELEMENT::STACK) {
                    const auto indices = feature->elements[ii].indices_;
                    auto last_pair = feature->sequence[indices[0]] + feature->sequence[indices[1]];

                    for(auto jj = 2; jj < indices.size(); ii += 2) {
                        const auto current_pair = feature->sequence[indices[jj]] + feature->sequence[indices[jj+1]];

                        if(last_pair == current_pair) {
                            score -= 1.f;
                        }
                        last_pair = current_pair;
                    }
                }
            }

            return score;
        }
    };

}
#endif // RNAMAKE_ELITWISTEDBASEPAIRS_H 