//
// Created by cjurich on 10/19/20.
//

#ifndef RNAMAKE_PCT_GC_PAIRS_H
#define RNAMAKE_PCT_GC_PAIRS_H

#include <base/types.h>

#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class PctGCPairs : public Strategy2D {


    public:
        PctGCPairs() : Strategy2D() {
            params_ = {88.50605747, 6.09781405};
            name_ = "PctGCPairs";
        }

    public:
        float
        score(Feature2DOP const & feature) {

            const auto pct_gc = feature->gc / (feature->gu + feature->gc + feature->ua);

            return params_[0] + pct_gc*params_[1];

        }
    };
}


#endif //RNAMAKE_PCT_GC_PAIRS_H
