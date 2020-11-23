#ifndef RNAMAKE_ELIENERGYLIMITINTETRALOOPS_H
#define RNAMAKE_ELIENERGYLIMITINTETRALOOPS_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class EliEnergyLimitInTetraloops : public Strategy2D {
    private:

    public:
        EliEnergyLimitInTetraloops() : Strategy2D() {
            params_ = {4.5, -10.};
            name_ = "eli_energy_limit_in_tetraloops";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);
            auto tetraloop_found(false);
            const auto& elements = feature->elements;

            for(const auto& element : elements) {
                if(element.type_ != RNAELEMENT::LOOP) {
                    continue;
                }
                const auto closing_pairs = element.get_loop_closing_pairs(
                        feature->sequence, feature->e_pairmap
                        );

                if(closing_pairs.size() == 1 && element.indices_.size() == 4) {
                    tetraloop_found = true;
                    if(element.score_ > params_[0]) {
                        result += params_[1];
                    }
                }
            }

            if(tetraloop_found) {
                return result;
            } else {
                return rnamake2d::UNSCORABLE;
            }
        }
    };

}
#endif // RNAMAKE_ELIENERGYLIMITINTETRALOOPS_H 