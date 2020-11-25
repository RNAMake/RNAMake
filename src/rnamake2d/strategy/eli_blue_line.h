#ifndef RNAMAKE_ELIBLUELINE_H
#define RNAMAKE_ELIBLUELINE_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class EliBlueLine : public Strategy2D {
    private:

    public:
        EliBlueLine() : Strategy2D() {
            params_ = { -1.f };
            name_ = "eli_blue_line";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);

            for(const auto& elem : feature->elements) {
                if (elem.type_ == rnamake2d::RNAELEMENT::STACK) {
                    auto first_string = String();
                    auto second_string = String();

                    for(auto jj = 0; jj < elem.indices_.size(); jj += 2) {
                        first_string += feature->sequence[elem.indices_[jj]];
                        second_string += feature->sequence[elem.indices_[jj + 1]];
                    }

                    std::transform(first_string.begin(), first_string.end(), first_string.begin(), toupper);
                    std::transform(second_string.begin(), second_string.end(), second_string.begin(), toupper);


                    for(auto jj = 0 ; jj < first_string.size(); ++jj) {
                        if(first_string.substr(jj).size() > 3 && first_string.substr(jj,4) == "UUUU"){
                            result += params_[0];
                        }
                        if(second_string.substr(jj).size() > 3 && second_string.substr(jj,4) == "UUUU"){
                            result += params_[0];
                        }
                    }
                }
            }
            return result;
        }
    };

}
#endif // RNAMAKE_ELIBLUELINE_H 