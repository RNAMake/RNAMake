#ifndef RNAMAKE_KKOHLITESTBYKKOHLI_H
#define RNAMAKE_KKOHLITESTBYKKOHLI_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class KkohliTestByKkohli : public Strategy2D {
    private:

    public:
        KkohliTestByKkohli() : Strategy2D() {
            params_ = {-0.408956766363662,
                       2.4557888979860456,
                       -4.4216292075931785,
                       -0.24488055252619767,
                       4.48164924686959,
                       40.31519149319874};
            name_ = "kkohli_test_by_kkohli";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);
            const auto& sequence = feature->sequence;
            auto g_count(0);

            for(auto ii = 0; ii < feature->sequence.size(); ++ii) {
                if(sequence[ii] == 'G')  {
                    ++g_count;
                }
            }
            const auto g_rate = float(g_count) / float(sequence.size());

            auto closing_pairs_count(0);
            auto closing_gc_count(0);

            for(auto ii = 0; ii < feature->elements.size(); ++ii) {
                if(feature->elements[ii].type_ == RNAELEMENT::LOOP) {
                    if(feature->elements[ii].score_ > params_[0]) {
                        result -= params_[1];
                    }

                    const auto closing_pairs = feature->elements[ii].get_loop_closing_pairs(
                            feature->sequence, feature->e_pairmap);
                    for(auto jj = 0; jj < closing_pairs.size(); ++jj){
                        ++closing_pairs_count;
                        if(closing_pairs[jj] == "GC" || closing_pairs[jj] == "CG") {
                            ++closing_gc_count;
                        }
                    }
                }
            }
            const auto closing_gc_rate = float(closing_gc_count) / float(closing_pairs_count);

            if(feature->meltpoint < 97.f || feature->meltpoint > 107.f)  {
                result -= params_[2];
            }

            if(feature->fe > -30.f)  {
                result -= (feature->fe - (-30.f ))*params_[3];
            }

            if(g_rate < 0.25f || g_rate > 0.35f)  {
                result -= params_[4];
            }

            if(closing_gc_rate < 0.6f) {
                result -= (0.6f - closing_gc_rate)*params_[5];
            }

            return result;
        }
    };

}
#endif // RNAMAKE_KKOHLITESTBYKKOHLI_H 