#ifndef RNAMAKE_NON_GC_CLOSING_H
#define RNAMAKE_NON_GC_CLOSING_H

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class JunctionNonGCClosing : public Strategy2D {
    public:
        JunctionNonGCClosing() : Strategy2D() {
            name_ = "JunctionNonGCClosing";
            params_ = {17.29684124};
        }

    public:
        float
        score(Feature2DOP const& feature) override {
            auto result(100.f);

            for(const auto& jnc : feature->junctions()) {
                result -= params_[0]*(jnc->au + jnc->gu + jnc->unknown);
            }

            return result;
        }
    };
}

#endif // RNAMAKE_NON_GC_CLOSING_H