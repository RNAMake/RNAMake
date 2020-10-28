#ifndef RNAMAKE_CJMISMATCH_H
#define RNAMAKE_CJMISMATCH_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class CJMismatch : public Strategy2D {
    private:
        unsigned long shift_limit_ = 3;
    public:
        CJMismatch() : Strategy2D() {
            params_ = { 7.55914918,  9.94522105, 12.58366418};
            name_ = "CJMismatch";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);

            for(auto shift = 0 ; shift < std::min(shift_limit_, feature->sequence.size()) ; ++shift) {
                for(auto index = 0; index < feature->sequence.size() - shift ; ++index) {
                    const auto pair = String{feature->sequence[index]}
                    + String{feature->sequence[feature->sequence.size() - index - shift - 1]};

                    if(pair == "AU" || pair == "UA") {
                        result -= params_[0];
                    } else if (pair == "UG" || pair == "GU") {
                        result -= params_[1];
                    } else if (pair == "GC" || pair == "CG") {
                        result -= params_[2];
                    }
                }
            }

            return result;
        }
    };

}
#endif // RNAMAKE_CJMISMATCH_H 