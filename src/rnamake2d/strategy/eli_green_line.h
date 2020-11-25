#ifndef RNAMAKE_ELIGREENLINE_H
#define RNAMAKE_ELIGREENLINE_H

#include <unordered_set>
#include <regex>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class EliGreenLine : public Strategy2D {
    private:

    public:
        EliGreenLine() : Strategy2D() {
            params_ = {-13.2821875};
            name_ = "eli_green_line";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);

            for(auto ii = 0; ii < feature->elements.size(); ++ii) {
                const auto& elem = feature->elements[ii];

                if(elem.type_ == RNAELEMENT::STACK) {
                    String first_string("");
                    String second_string("");

                    for(auto jj = 0; jj < elem.indices_.size(); jj += 2) {
                        first_string += feature->sequence[elem.indices_[jj]];
                        second_string += feature->sequence[elem.indices_[jj + 1]];
                    }

                    std::transform(first_string.begin(), first_string.end(), first_string.begin(), toupper);
                    std::transform(second_string.begin(), second_string.end(), second_string.begin(), toupper);

                    for(auto jj = 0; jj < first_string.size(); ++jj) {
                        if(first_string.substr(jj).size() > 2 && first_string.substr(jj, 3) == "CCC") {
                            result += params_[0];
                        }
                        if(second_string.substr(jj).size() > 2 && second_string.substr(jj, 3) == "CCC") {
                            result += params_[0];
                        }
                    }
                }
            }
            return result;
        }
    };

}
#endif // RNAMAKE_ELIGREENLINE_H 