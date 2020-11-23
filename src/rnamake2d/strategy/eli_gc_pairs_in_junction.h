#ifndef RNAMAKE_ELIGCPAIRSINJUNCTION_H
#define RNAMAKE_ELIGCPAIRSINJUNCTION_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class EliGCPairsInJunction : public Strategy2D {
    private:

    public:
        EliGCPairsInJunction() : Strategy2D() {
            params_ = {1.f};
            name_ = "eli_gc_pairs_in_junction";
        }

    private:
        static bool
        is_neighbor_paried( const int idx, Ints const& pairmap, const int length) {
            if(((idx > 0) && (pairmap[idx - 1] < 0)) or (
                    (idx < length - 1)  && (pairmap[idx + 1] < 0 )
                    ))  {
                return true;
            } else {
                return false;
            }
        }

    private:
        static bool
        is_gc(const char base_a, const char base_b)  {
            if( (base_a == 'G' && base_b == 'C') || (base_a == 'C' && base_b == 'G') ) {
                return true;
            } else {
                return false;
            }
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);

            auto junct_ct(0);
            auto junct_gc(0);

            for(auto ii = 0; ii < feature->e_pairmap.size(); ++ii) {
                const auto idx = feature->e_pairmap[ii];

                if(idx > ii) {
                    if(!is_neighbor_paried(ii, feature->e_pairmap, feature->e_pairmap.size()) ||
                            !is_neighbor_paried(idx, feature->e_pairmap, feature->e_pairmap.size()) ) {
                        ++junct_ct;
                        if(is_gc(feature->sequence[ii], feature->sequence[idx])) {
                            ++junct_gc;
                        }
                    }

                    result += -params_[0]*float(junct_ct - junct_gc);
                }

            }
            return result;
        }
    };

}
#endif // RNAMAKE_ELIGCPAIRSINJUNCTION_H 