#ifndef RNAMAKE_ELIREDLINE_H
#define RNAMAKE_ELIREDLINE_H

#include <algorithm>
#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class EliRedLine : public Strategy2D {
    private:

    public:
        EliRedLine() : Strategy2D() {
            params_ = {-9.748725585937487};
            name_ = "eli_red_line";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);

            for(auto ii = 0; ii < feature->elements.size(); ++ii) {
                const auto& elem = feature->elements[ii];

                if(elem.type_ == RNAELEMENT::STACK) {
                   auto first_string = String{};
                   auto second_string = String{};

                   for(auto jj = 0; jj < elem.indices_.size(); jj += 2 ) {
                         first_string += feature->sequence[elem.indices_[jj]];
                         second_string += feature->sequence[elem.indices_[jj + 1]];
                    }
                    const auto banned = String{"GGG"};

                    std::transform(first_string.begin(), first_string.end(), first_string.begin(), toupper);
                    std::transform(second_string.begin(), second_string.end(), second_string.begin(), toupper);

                    for(auto jj = 0; jj < first_string.size(); ++jj) {
                        if(first_string.substr(jj).find(banned) != String::npos) {
                            result += params_[0];
                        }
                        if(second_string.substr(jj).find(banned) != String::npos) {
                            result += params_[0];
                        }
                    }
                }
            }
            return result;
        }
    };

}
#endif // RNAMAKE_ELIREDLINE_H 