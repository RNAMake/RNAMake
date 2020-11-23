#ifndef RNAMAKE_RNAMAKEALDOMISMATCH_H
#define RNAMAKE_RNAMAKEALDOMISMATCH_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class RNAMakeAldoMismatch : public Strategy2D {
    private:

    public:
        RNAMakeAldoMismatch() : Strategy2D() {
            params_ = {1, 1, 1};
            name_ = "AldoMismatch";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);

            const auto length = feature->target.size();

            const auto limit = length < SHIFT_LIMIT ? length : SHIFT_LIMIT;

            for(auto shift = 0; shift < limit ; ++shift) {
                for(auto index = 0; index < (length - shift) ; ++index )  {
                    auto pair = String(2, ' ');
                    pair[0] = feature->sequence[index];
                    pair[1] = feature->sequence[length - index - shift - 1];

                    if(pair == "AU" || pair == "UA") {
                        result -= params_[0]/float(length);
                    } else if (pair == "UG" || pair == "GU") {
                        result -= params_[1]/float(length);
                    } else if (pair == "GC" || pair == "CG") {
                        result -= params_[2]/float(length);
                    }
                }
            }
            return result;
        }
    };

}
#endif // RNAMAKE_RNAMAKEALDOMISMATCH_H 