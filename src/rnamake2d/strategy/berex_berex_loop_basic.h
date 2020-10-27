#ifndef RNAMAKE_BEREXBEREXLOOPBASIC_H
#define RNAMAKE_BEREXBEREXLOOPBASIC_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class BerexBerexLoopBasic : public Strategy2D {
    private:

    public:
        BerexBerexLoopBasic() : Strategy2D() {
            params_ = {4.8526841782844325,
                       3.0696390175733637,
                       1.9282212756410848,
                       -1.1981440815216913,
                       16.70193159320271,
                       -18.661384058745625,
                       0.99512254095277,
                       -7.899441097535717
            };
            name_ = "berex_berex_loop_basic";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(80.f);

            for(auto ii = 0; ii < feature->elements.size(); ++ii) {
                const auto elem = feature->elements[ii];

                if(elem.type_ == RNAELEMENT::LOOP) {
                    continue;
                }

                const auto loop_groups = elem.get_loop_groups();
                const auto closing_pairs = elem.get_loop_closing_pairs(feature->sequence, feature->pairmap);

                const auto indices = elem.indices_;

                if(loop_groups.size() == 1 && closing_pairs.size() == 2) {
                    if(closing_pairs[0] == "GC"
                        || closing_pairs[0] == "CG"
                        || closing_pairs[1] == "GC"
                        || closing_pairs[1] == "CG") {
                        result += params_[0];
                    }
                } else if ( loop_groups.size() == 2 && closing_pairs.size() == 2) {
                    if(loop_groups[0].size() == 1 && loop_groups[1].size() ==1 ) {
                       if(feature->sequence[loop_groups[0][0]] == 'G'
                            && feature->sequence[loop_groups[1][0]] == 'G') {
                           result += params_[1];
                       }
                    } else if (loop_groups[0].size() + loop_groups[1].size() == 3) {
                        auto bottom(0), top(0);
                        if(loop_groups[0].size() == 1) {
                            bottom = 0;
                            top = 1;
                        } else {
                            bottom = 1;
                            top = 0;
                        }
                        if(feature->sequence[loop_groups[bottom][0]] == 'G') {
                            result += params_[2];
                        }
                    } else if (loop_groups[0].size() + loop_groups[1].size() == 4) {
                        if(loop_groups[0].size() == 2 && loop_groups[1].size() == 2) {
                            if ( (feature->sequence[loop_groups[0][0]] == 'G'
                                 or feature->sequence[loop_groups[0][1]] == 'G')
                                 and
                                 (feature->sequence[loop_groups[1][0]] == 'G'
                                 or feature->sequence[loop_groups[1][1]] == 'G')
                            ) {
                                result += params_[3];
                            }
                        } else {
                            auto bottom(0), top(0);
                            if(loop_groups[0].size() == 1) {
                                bottom = 0;
                                top = 1;
                            } else {
                                bottom = 1;
                                top = 0;
                            }
                            auto top_g_count(0);

                            if(feature->sequence[loop_groups[top][0]]  == 'G') {
                                ++top_g_count;
                            }
                            if(feature->sequence[loop_groups[top][1]]  == 'G') {
                                ++top_g_count;
                            }
                            if(feature->sequence[loop_groups[top][2]]  == 'G') {
                                ++top_g_count;
                            }
                            if(feature->sequence[loop_groups[bottom][0]]  == 'G') {
                                result += params_[4];
                            }
                            if(top_g_count == 2) {
                                result += params_[4];
                            }
                        }
                    } else {
                        if(loop_groups.size() == 1
                        && loop_groups[0].size() == 4
                        && closing_pairs.size() == 1) {
                            if(feature->sequence[loop_groups[0][0]] == 'G') {
                                result += params_[5];
                            }
                        } else {
                            auto g_count(0);
                            for(auto jj = 0; jj < indices.size(); ++jj) {
                                if(feature->sequence[indices[jj]] == 'G') {
                                    ++g_count;
                                }
                            }

                            if(g_count == 1) {
                                result += params_[6];
                            } else if (g_count >= 2) {
                                result += params_[7];
                            }
                        }
                    }
                }
            }

            return result;
        }
    };

}
#endif // RNAMAKE_BEREXBEREXLOOPBASIC_H 