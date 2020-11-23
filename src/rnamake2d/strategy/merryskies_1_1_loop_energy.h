#ifndef RNAMAKE_MERRYSKIES_1_1_LOOPENERGY_H
#define RNAMAKE_MERRYSKIES_1_1_LOOPENERGY_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class MerrySkies_1_1_LoopEnergy : public Strategy2D {
    private:

    public:
        MerrySkies_1_1_LoopEnergy() : Strategy2D() {
            params_ = {-0.39, 5.f};
            name_ = "merryskies_1_1_loop_energy";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            const auto& elements = feature->elements;
            const auto& sequence = feature->sequence;
            const auto& pairmap = feature->e_pairmap;

            auto score(80.f);
            auto count(0);
            for(auto ii = 0; ii < elements.size(); ++ii) {
                const auto& elem = elements[ii];
                if(elem.type_ != RNAELEMENT::LOOP) {
                    continue;
                }

                const auto loop_groups = elem.get_loop_groups();
                const auto closeing_pairs = elem.get_loop_closing_pairs(sequence, pairmap);

                if(loop_groups.size() == 2 && closeing_pairs.size() == 2) {
                    if(loop_groups[0].size() == 1 && loop_groups[1].size() == 1) {
                        count += 1;
                        if(elements[ii].score_ < params_[0]) {
                            score += params_[1];
                        }
                    }
                }
            }

            if(count > 0)  {
                return score;
            } else {
                return UNSCORABLE;
            }
        }
    };

}
#endif // RNAMAKE_MERRYSKIES_1_1_LOOPENERGY_H 