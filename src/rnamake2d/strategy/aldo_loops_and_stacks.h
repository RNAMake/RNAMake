#ifndef RNAMAKE_ALDOLOOPSANDSTACKS_H
#define RNAMAKE_ALDOLOOPSANDSTACKS_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class AldoLoopsAndStacks : public Strategy2D {
    private:

    public:
        AldoLoopsAndStacks() : Strategy2D() {
            params_ = {0.21599778719663132,
                       0.2461934637984866,
                       1.6804238484498124,
                       1.3549843846144776,
                       0.9987290691132436,
                       0.6308787845195161,
                       0.08692816268745038,
                       0.6606132364847958,
                       273.3865485603875};
            name_ = "aldo_loops_and_stacks";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(0.f);

            for(auto ii = 0; ii < feature->elements.size(); ++ii) {
                const auto& elem = feature->elements[ii];
                if(elem.type_ == RNAELEMENT::STACK) {
                    const auto& stack_len =  elem.get_stack_length();
                    for(auto jj = 0; jj < stack_len; ++jj) {
                        const auto pair = elem.get_pair_from_stack(jj, feature->sequence);
                        if(pair == "GC" || pair == "CG") {
                            if(jj == 0 || jj == stack_len - 1) {
                                result += params_[0];
                            } else {
                                result += params_[1];
                            }
                        } else if (pair == "UA" || pair == "AU") {
                            if(jj == 0 || jj == stack_len - 1) {
                                result += params_[2];
                            } else {
                                result += params_[3];
                            }
                        } else {
                            result += params_[4];
                        }
                    }
                }
            }

            for(auto ii = 0 ; ii < feature->sequence.size(); ++ii) {
                if(feature->target[ii] == '.' && feature->sequence[ii] == 'A') {
                    result += params_[5];
                } else if (feature->target[ii] == '.' && feature->sequence[ii] == 'G') {
                    result += params_[6];
                }
            }

            const auto modifier = 1.f - std::abs(feature->gc/(feature->gc + feature->ua + feature->gu) - params_[7]);
            return params_[8]*modifier*result / float(feature->sequence.size());
        }
    };

}
#endif // RNAMAKE_ALDOLOOPSANDSTACKS_H 