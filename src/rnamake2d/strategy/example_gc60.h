#ifndef RNAMAKE_EXAMPLEGC60_H
#define RNAMAKE_EXAMPLEGC60_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class ExampleGC60 : public Strategy2D {
    private:

    public:
        ExampleGC60() : Strategy2D() {
            params_ = {0.5746352738846094};
            name_ = "example_gc60";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            return 100.f - std::abs(params_[0] - feature->gc / (feature->gc + feature->gu + feature->ua))*100.f;
        }
    };

}
#endif // RNAMAKE_EXAMPLEGC60_H 