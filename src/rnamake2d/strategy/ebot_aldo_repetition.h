#ifndef RNAMAKE_EBOTALDOREPETITION_H
#define RNAMAKE_EBOTALDOREPETITION_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class EbotAldoRepetition : public Strategy2D {
    private:

    public:
        EbotAldoRepetition() : Strategy2D() {
            params_ = {0, 0, 5, 5, 5, 5};
            name_ = "aldo_repetition";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            const auto& sequence = feature->sequence;
            const auto length = sequence.size();
            auto rep_cnt = Ints(6);
            auto n_gram = std::map<String,int>();

            for(auto ii = 0; ii < length; ++ii)  {
                for(auto jj = 2; jj < 6; ++jj) {
                    if( ii + jj >= length) {
                        break;
                    }
                    const auto substring = sequence.substr(ii, jj);

                    if(n_gram.find(substring) != n_gram.end()) {
                        n_gram[substring]++;
                    } else {
                        n_gram[substring] = 0;
                    }
                }
            }

            auto result(100.f);

            for(auto& [substring, ct] : n_gram)  {
                rep_cnt[substring.size()] += ct;
            }
            for(auto ii = 1; ii < 6; ++ii) {
                result -= float(params_[ii] * rep_cnt[ii]) / length;
            }

            return result;
        }
    };

}
#endif // RNAMAKE_EBOTALDOREPETITION_H 