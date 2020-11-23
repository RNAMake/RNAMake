#ifndef RNAMAKE_DINGTETRALOOPPATTERN_H
#define RNAMAKE_DINGTETRALOOPPATTERN_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class DingTetraloopPattern : public Strategy2D {
    private:

    public:
        DingTetraloopPattern() : Strategy2D() {
            params_ = {-5.};
            name_ = "ding_tetraloop_pattern";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);

            const auto& sequence = feature->sequence;
            const auto& elements = feature->elements;

            for(const auto& element : elements)  {
                if(element.type_ != RNAELEMENT::LOOP) {
                    continue;
                }

                if(element.get_loop_closing_pairs(feature->sequence, feature->e_pairmap).size() == 1
                    and
                    element.indices_.size() == 4
                )  {

                    auto tetraloop_string = String(4, ' ');
                    tetraloop_string[0] = sequence[element.indices_[0]];
                    tetraloop_string[1] = sequence[element.indices_[1]];
                    tetraloop_string[2] = sequence[element.indices_[2]];
                    tetraloop_string[3] = sequence[element.indices_[3]];

                    if(tetraloop_string != "AAAA")  {
                        result += params_[0];
                    }
                }
            }
            return result;
        }
    };

}
#endif // RNAMAKE_DINGTETRALOOPPATTERN_H 