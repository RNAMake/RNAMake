#ifndef RNAMAKE_PENGUIANCLEANDOTPLOT_H
#define RNAMAKE_PENGUIANCLEANDOTPLOT_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class PenguianCleanDotplot : public Strategy2D {
    private:

    public:
        PenguianCleanDotplot() : Strategy2D() {
            params_ = {};
            name_ = "penguian_clean_dotplot";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto penalty(0.f);
            const auto npairs = feature->gu + feature->gc + feature->ua;
            auto& dotplot = feature->dotplot;
            for(auto ii = 0; ii < feature->dotplot.size(); ++ii) {
                const auto i_index = dotplot[ii].i;
                const auto j_index = dotplot[ii].j;

                if(feature->pairmap[i_index] != j_index) {
                    penalty += dotplot[ii].p;
                }
            }
            return 100.f - (penalty/npairs);
        }
    };

}
#endif // RNAMAKE_PENGUIANCLEANDOTPLOT_H 