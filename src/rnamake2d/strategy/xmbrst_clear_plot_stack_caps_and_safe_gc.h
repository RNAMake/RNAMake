#ifndef RNAMAKE_XMBRSTCLEARPLOTSTACKCAPSANDSAFEGC_H
#define RNAMAKE_XMBRSTCLEARPLOTSTACKCAPSANDSAFEGC_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class XmbrstClearPlotStackCapsAndSafeGC : public Strategy2D {
    private:

    public:
        XmbrstClearPlotStackCapsAndSafeGC() : Strategy2D() {
            params_ = {1.0, 1.0, 2.0, 0.8};
            name_ = "xmbrst_clear_plot_stack_caps_and_safe_gc";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto penalty(0.f);
            const auto n = feature->sequence.size();
            const auto& sequence = feature->sequence;
            const auto npairs = feature->gc + feature->gu + feature->ua;
            const auto& dotplot = feature->dotplot;
            const auto& pairmap = feature->e_pairmap;
            const auto& elements = feature->elements;

            for(auto ii = 0; ii < dotplot.size(); ++ii) {
                const auto i_index = dotplot[ii].i;
                const auto j_index = dotplot[ii].j;

                if(pairmap[i_index] != j_index) {
                    penalty += dotplot[ii].p;
                }
            }
            auto plotscore(0.f);
            if(npairs > 0) {
                plotscore = 1.0f - penalty/npairs;
            } else {
                plotscore = 0.f;
            }
            auto gc_penalty(0);
            if (npairs > 0.f) {
                if(feature->gc/npairs > params_[3]) {
                    gc_penalty = 1;
                }
            }

            auto cap_score(0.f);
            auto stack_count(0);
            for(auto ii = 0; ii < elements.size(); ++ii) {
                if(elements[ii].type_ == RNAELEMENT::STACK) {
                    const auto stacklen = int(elements[ii].indices_.size())/2;
                    stack_count++;
                    if(stacklen == 1) {
                        const auto pair = elements[ii].get_pair_from_stack(0, sequence);
                        if(pair == "GC" || pair == "CG") {
                            cap_score++;
                        }
                    } else if (stacklen == 2) {
                        auto pair = elements[ii].get_pair_from_stack(0, sequence);
                        if(pair == "GC" || pair == "CG") {
                            cap_score += 0.5;
                        }
                        pair = elements[ii].get_pair_from_stack(1, sequence);
                        if(pair == "CG" || pair == "GC") {
                            cap_score += 0.5;
                        }
                    } else if (stacklen == 3) {
                        for(auto ii = 0; ii < 3; ++ii)  {
                            const auto pair = elements[ii].get_pair_from_stack(ii, sequence);
                            if(pair == "GC" || pair == "CG") {
                                cap_score += 0.4;
                            }
                        }
                    } else {
                        auto pair = elements[ii].get_pair_from_stack(0, sequence);
                        if(pair == "GC" || pair == "CG") {
                            cap_score += 1.f / 3.0f;
                        }
                        pair = elements[ii].get_pair_from_stack(0, sequence);
                        if(pair == "GC" || pair == "CG") {
                            cap_score += 1.f / 6.0f;
                        }
                        pair = elements[ii].get_pair_from_stack(stacklen - 2, sequence);
                        if(pair == "GC" || pair == "CG") {
                            cap_score += 1.f / 6.0f;
                        }
                        pair = elements[ii].get_pair_from_stack(stacklen - 1, sequence);
                        if(pair == "GC" || pair == "CG") {
                            cap_score += 1.f / 3.0f;
                        }

                    }
                }
            }

            if (stack_count > 0) {
                cap_score = float(cap_score) / float(stack_count);
            }

            return (2.0f + cap_score*params_[1] + plotscore*params_[0] - gc_penalty*params_[2])*25.f;
        }
    };

}
#endif // RNAMAKE_XMBRSTCLEARPLOTSTACKCAPSANDSAFEGC_H 