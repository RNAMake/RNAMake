#ifndef RNAMAKE_RNAMAKEALDOLOOPSANDSTACKS_H
#define RNAMAKE_RNAMAKEALDOLOOPSANDSTACKS_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class RNAMakeAldoLoopsAndStacks : public Strategy2D {
    private:

    public:
        RNAMakeAldoLoopsAndStacks() : Strategy2D() {
            params_ = {1, 1, 1, 1, 1, 1, 1, 1, 1};
            name_ = "AldoLoopAndStacks";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto gc_non_closing = 0;
            auto gc_closing_pair = 0;
            auto ua_non_closing = 0;
            auto ua_closing_pair = 0;
            auto other_pairs = 0;
            auto a_unpaired = 0;
            auto g_unpaired = 0;
            auto gc_pairs = 0;
            auto number_pairs = 0;
            const auto structure_length = feature->target.size();

            for(const auto& stack : feature->elements)  {
                if(stack.type_ != RNAELEMENT::STACK) {
                    continue;
                }
                auto [lhs, rhs] = stack.get_sides(feature->sequence);
                const auto sequence = lhs + "&" + rhs;
                for(auto index = 0; index < stack.get_stack_length(); ++index) {
                    number_pairs += 1;
                    auto pair = String(2, ' ');
                    pair[0] = sequence[index];
                    pair[1] = *(sequence.rbegin()+ index);
                    assert(pair.find('&') == String::npos);

                    if( pair == "GC" or pair == "CG") {
                        gc_pairs += 1;
                        if (index == 0 or index == stack.get_stack_length() - 1) {
                            gc_closing_pair += 1;
                        } else {
                            gc_non_closing += 1;
                        }
                    } else if (pair == "AU" or pair == "UA") {
                        if ( index == 0 or index == stack.get_stack_length() - 1) {
                            ua_closing_pair += 1;
                        } else {
                            ua_non_closing += 1;
                        }
                    } else {
                        other_pairs += 1;
                    }

                }

            }

            for(auto ii = 0; ii < feature->sequence.size(); ++ii) {
                const auto db = feature->target[ii];
                const auto nt = feature->sequence[ii];
                if(db == '.') {
                    if(nt == 'A') {
                        a_unpaired += 1  ;
                    } else if ( nt == 'G' ){
                        g_unpaired += 1;
                    }
                }
            }

            number_pairs = std::count(feature->target.begin(), feature->target.end(), '(');
            const auto modifier = 1.0f - abs(float(gc_pairs/number_pairs) - params_[7] );

            const auto score = params_[0]*gc_closing_pair + params_[1]*gc_non_closing + params_[2]*ua_non_closing
                    + params_[3]*ua_non_closing + params_[4]*other_pairs + params_[5]*a_unpaired + params_[6]*g_unpaired;

            return params_[8]*modifier*score/float(structure_length);
        }
    };

}
#endif // RNAMAKE_RNAMAKEALDOLOOPSANDSTACKS_H 