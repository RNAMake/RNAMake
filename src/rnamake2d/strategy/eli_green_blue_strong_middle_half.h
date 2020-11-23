#ifndef RNAMAKE_ELIGREENBLUESTRONGMIDDLEHALF_H
#define RNAMAKE_ELIGREENBLUESTRONGMIDDLEHALF_H

#include <regex>
#include <algorithm>
#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class EliGreenBlueStrongMiddleHalf : public Strategy2D {
    private:

    public:
        EliGreenBlueStrongMiddleHalf() : Strategy2D() {
            params_ = {-0.5749999999999797};
            name_ = "eli_green_blue_strong_middle_half";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto score(100.f);
            for(auto ii = 0; ii < feature->elements.size(); ++ii) {
                const auto& elem = feature->elements[ii];

                if(elem.type_ == RNAELEMENT::STACK) {
                    String first_string("");
                    String second_string("");

                    for(auto jj = 0 ; jj < elem.indices_.size(); jj += 2) {
                        first_string += feature->sequence[elem.indices_[jj]];
                        second_string += feature->sequence[elem.indices_[jj+1]];
                    }

                    std::transform(first_string.begin(), first_string.end(), first_string.begin(), toupper);
                    std::transform(second_string.begin(), second_string.end(), second_string.begin(), toupper);

                    auto banned = std::regex();
                    if(elem.indices_.size() < 12) {
                        banned = std::regex("^[U|C][U|C][U|C][U|C]");
                    } else {
                        banned = std::regex("^[U|C][U|C][U|C][U|C][U|C]");
                    }

                    for(auto jj = 0; jj < first_string.size(); ++jj) {
                        if(!findall(first_string.substr(jj),banned).empty()) {
                            score += params_[0];
                        }
                        if(!findall(second_string.substr(jj),banned).empty()) {
                            score += params_[0];
                        }
                    }
                }
            }
            return score;
        }
    };

}
#endif // RNAMAKE_ELIGREENBLUESTRONGMIDDLEHALF_H 