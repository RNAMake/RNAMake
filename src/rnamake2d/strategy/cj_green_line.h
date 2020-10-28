#ifndef RNAMAKE_CJGREENLINE_H
#define RNAMAKE_CJGREENLINE_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class CJGreenLine : public Strategy2D {
    private:

    public:
        CJGreenLine() : Strategy2D() {
            params_ = {6.84687915};
            name_ = "CJGreenLine";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);

            for(auto start = 0; start < feature->sequence.size() - 3 ; ++start) {
                const auto structure = feature->structure.substr( start, 3);
                if(structure.find('.') != String::npos) {
                    continue;
                }

                const auto sequence = feature->sequence.substr( start, 3);
                if(sequence == "CCC") {
                    result -= params_[0];
                }
            }

            return result;
        }
    };

}
#endif // RNAMAKE_CJGREENLINE_H 