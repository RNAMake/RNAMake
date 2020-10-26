#ifndef RNAMAKE_CLOLLINGSINPLACE_H
#define RNAMAKE_CLOLLINGSINPLACE_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class ClollinGsInPlace : public Strategy2D {
    private:

    public:
        ClollinGsInPlace() : Strategy2D() {
            params_ = {};
            name_ = "";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(80.f);
            auto count(0);
            for(auto ii = 0 ; ii < feature->elements.size(); ++ii) {
                if(feature->elements[ii].type_ == RNAELEMENT::LOOP) {
                    auto loop_groups = feature->elements[ii].get_loop_groups();
                    if(loop_groups.size() == 1 && loop_groups[0].size() > 3) {
                       ++count;
                       if(feature->sequence[loop_groups[0][0]] == 'G') {
                            result += 1.f;
                       }
                    }
                }
            }

            if(count < 1) {
                return UNSCORABLE;
            }

            return result;
        }
    };

}
#endif // RNAMAKE_CLOLLINGSINPLACE_H 