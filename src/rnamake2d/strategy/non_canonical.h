#ifndef RNAMAKE_NON_CANONICAL_H
#define RNAMAKE_NON_CANONICAL_H

#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class NonCanonical : public Strategy2D {
    public:
        NonCanonical() : Strategy2D() {
            name_ = "NonCanonical";
            params_ = {1.20440825};
        }

    public:
        float
        score(Feature2DOP const & feature ) override {
            const auto pairs = std::count(feature->structure.begin(), feature->structure.end(),')');
            return 100.f - (pairs - feature->gu - feature->gc - feature->ua)*params_[0];
        }
    };
}



#endif // RNAMAKE_NON_CANONICAL_H
