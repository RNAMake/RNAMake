#ifndef RNAMAKE_MAT747_31_LOOPS_H
#define RNAMAKE_MAT747_31_LOOPS_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class Mat747_31_Loops : public Strategy2D {
    private:

    public:
        Mat747_31_Loops() : Strategy2D() {
            params_ = {1};
            name_ = "mat747_31_loops";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            const auto& elements = feature->elements;
            auto loop_count(0), right_loop_count(0);
            const auto& sequence = feature->sequence;
            const auto& pairmap = feature->e_pairmap;

            for(auto ii = 0; ii < elements.size(); ++ii) {
                const auto& element = elements[ii];

                if(element.type_ == RNAELEMENT::LOOP) {
                    const auto& loop_groups = elements[ii].get_loop_groups();
                    if(loop_groups.size() == 2) {
                        const auto grp_1_len = loop_groups[0].size();
                        const auto grp_2_len = loop_groups[1].size();
                        const auto closeing_pairs = element.get_loop_closing_pairs(sequence, pairmap);

                        if(grp_1_len == 1 && grp_2_len == 3)  {
                            loop_count += 0;
                        } else if (grp_1_len == 3 && grp_2_len == 1) {
                            loop_count += 1;
                            if(closeing_pairs[0] == "GC" and *closeing_pairs.rbegin() == "GC") {
                                right_loop_count += 1;
                            }
                        }

                    }
                }
            }

            if(loop_count == 0) {
                return rnamake2d::UNSCORABLE;
            }

            return 100.f - float(loop_count - right_loop_count)*params_[0];
    }
    };

}
#endif // RNAMAKE_MAT747_31_LOOPS_H 