#ifndef RNAMAKE_MERRYSKIES_1_1_LOOP_H
#define RNAMAKE_MERRYSKIES_1_1_LOOP_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class MerrySkies_1_1_Loop : public Strategy2D {
    private:

    public:
        MerrySkies_1_1_Loop() : Strategy2D() {
            params_ = {2.49909690821, -0.890563715395, 2.83384425139};
            name_ = "merryskies_1_1_loop";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(80.f);
            for(auto ii = 0; ii < feature->elements.size(); ++ii) {

                auto& elem = feature->elements[ii];
                if( elem.type_ != RNAELEMENT::LOOP) {
                    continue;
                }

                auto loop_groups = elem.get_loop_groups();
                auto closing_pairs = elem.get_loop_closing_pairs(feature->sequence, feature->pairmap);

                const auto& sequence = feature->sequence;
                if( loop_groups.size() == 2 and closing_pairs.size() == 2 ) { // case 2,3, 4
                    if (loop_groups[0].size() == 1 and loop_groups[1].size() == 1) {  // case
                        if (
                                sequence[loop_groups[0][0]] == 'G'
                                and sequence[loop_groups[1][0]] == 'G'
                                ) {
                            result += params_[0];
                        } else if (
                                sequence[loop_groups[0][0]] == 'G'
                                and sequence[loop_groups[1][0]] == 'A'
                        ) {
                            result += params_[1];
                        } else if (
                                sequence[loop_groups[0][0]] == 'A'
                                and sequence[loop_groups[1][0]] == 'G'
                        ) {
                            result += params_[1];
                        } else if (
                                sequence[loop_groups[0][0]] == 'A'
                                and sequence[loop_groups[1][0]] == 'A'
                        ) {
                            result+= params_[2];
                        }
                    }
                }
            }
            return result;
        }
    };

}
#endif // RNAMAKE_MERRYSKIES_1_1_LOOP_H 