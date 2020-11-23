#ifndef RNAMAKE_AREPEATS_H
#define RNAMAKE_AREPEATS_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class ARepeats : public Strategy2D {
    private:

    public:
        ARepeats() : Strategy2D() {
            params_ = {1};
            name_ = "ARepeats";
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
                auto forbidden = std::regex("^AAA");
                for(auto start = 0; start < stack.get_stack_length(); ++start) {
                    if (!findall(lhs.substr(start), forbidden).empty()) {
                        result -= params_[0];
                    }
                    if (!findall(rhs.substr(start), forbidden).empty()) {
                        result -= params_[0];
                    }
                }
            }

            return result;
        }
    };

}
#endif // RNAMAKE_AREPEATS_H 