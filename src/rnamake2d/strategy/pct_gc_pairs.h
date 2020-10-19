//
// Created by cjurich on 10/19/20.
//

#ifndef RNAMAKE_PCT_GC_PAIRS_H
#define RNAMAKE_PCT_GC_PAIRS_H

#include <base/types.h>

#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class PctGCPairs : public Strategy2D {
        float base_score = 88.50605747;
        float multiplier = 6.09781405;

    public:
        PctGCPairs() : Strategy2D() {
            name_ = "PctGCPairs";
        }

    public:
        float
        score(Feature2DOP const & feature) {

            const auto pct_gc = feature->gc / (feature->gu + feature->gc + feature->ua);

            return base_score + pct_gc*multiplier;

        }
    };
}


#endif //RNAMAKE_PCT_GC_PAIRS_H
