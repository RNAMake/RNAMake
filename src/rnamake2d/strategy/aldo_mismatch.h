#ifndef RNAMAKE_ALDOMISMATCH_H
#define RNAMAKE_ALDOMISMATCH_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>
#include <rnamake2d/rna_element.h>

namespace rnamake2d {
    class AldoMismatch : public Strategy2D {
    private:

    public:
        AldoMismatch() : Strategy2D() {
            params_ = {2.8295833320877066, 10.774596546216415, 10.98943877153096};
            name_ = "aldo_mismatch";
        }
    private:
        static
        int
        which_pair(char base_a, char base_b) {
            if(base_a > base_b) {
                std::swap(base_a, base_b);
            }
            if(base_a == 'A' && base_b == 'U') {
                return 0;
            } else if (base_a == 'G' && base_b == 'U') {
                return 1;
            } else if (base_a == 'C' && base_b == 'G') {
                return 2;
            } else {
                return -1;
            }
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            const int length = feature->pairmap.size();
            auto pair_cnt = Ints(PAIR_TYPES,0);
            auto result(100.f);
            for(auto shift = 0; shift < std::min(SHIFT_LIMIT, length); ++shift) {
                for(auto ii = 0; ii < (length - shift) ; ++ii) {
                    auto base_a = feature->sequence[ii];
                    auto base_b = feature->sequence[length - ii - shift - 1];
                    auto pair_num = which_pair(base_a, base_b);
                    if (pair_num >= 0) {
                        pair_cnt[pair_num] += 1;
                    }
                }

            }
            for(auto ii = 0; ii < PAIR_TYPES; ++ii ) {
                result -= (params_[ii] * pair_cnt[ii]) / length;
            }
            return result;
        }
    };

}
#endif // RNAMAKE_ALDOMISMATCH_H 