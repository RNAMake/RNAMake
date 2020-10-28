#ifndef RNAMAKE_ELINOBLUENUCLEOTIDESINHOOK_H
#define RNAMAKE_ELINOBLUENUCLEOTIDESINHOOK_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class EliNoBlueNucleotidesInHook : public Strategy2D {
    private:

    public:
        EliNoBlueNucleotidesInHook() : Strategy2D() {
            params_ = {-0.17043017578125408};
            name_ = "eli_no_blue_nucleotides_in_hook";
        }

    public:
        static
        int
        find_nonnegative_value(Ints const& array)  {
            auto index(0);
            while(array[index] == -1) {
                ++index;
            }

            return index;
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto first_pair_index = find_nonnegative_value(feature->e_pairmap);
            auto result(100.f);
            const auto last_pair_index = feature->pairmap[first_pair_index];

            for(auto ii = 0; ii < first_pair_index + 1; ++ii) {
                if(feature->sequence[ii] == 'U') {
                    result += params_[0];
                }
            }

            for(auto ii = last_pair_index; ii < feature->sequence.size() ; ++ii ) {
                if(feature->sequence[ii] == 'U') {
                    result += params_[0];
                }
            }

            return result;
        }
    };

}
#endif // RNAMAKE_ELINOBLUENUCLEOTIDESINHOOK_H 