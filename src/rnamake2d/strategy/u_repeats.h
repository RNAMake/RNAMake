#ifndef RNAMAKE_UREPEATS_H
#define RNAMAKE_UREPEATS_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class URepeats : public Strategy2D {
    private:

    public:
        URepeats() : Strategy2D() {
            params_ = {1.};
            name_ = "URepeats";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);

            for(const auto& stack : feature->elements) {
                if(stack.type_ != RNAELEMENT::STACK) {
                    continue;
                }
                auto [rhs, lhs] = stack.get_sides(feature->sequence);
                for(auto start = 0; start < stack.get_stack_length(); ++start) {
                    if (lhs.substr(start).size() > 2 && lhs.substr(start,3) == "UUU") {
                        result -= params_[0];
                    }
                    if (rhs.substr(start).size() > 2 && rhs.substr(start,3) == "UUU") {
                        result -= params_[0];
                    }
                }
            }

            return result;
        }
    };

}
#endif // RNAMAKE_UREPEATS_H 