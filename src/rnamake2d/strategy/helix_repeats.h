#ifndef RNAMAKE_HELIXREPEATS_H
#define RNAMAKE_HELIXREPEATS_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class HelixRepeats : public Strategy2D {
    private:

    public:
        HelixRepeats() : Strategy2D() {
            params_ = {};
            name_ = "HelixRepeats";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);

            for(auto& motif : feature->motifs) {
                if(motif->mtype() != util::MotifType::HELIX || motif->token() == "helix1") {
                    continue;
                }
                const auto& tk = motif->token();

                if(tk == "helix2") {
                   result -= params_[0]*num_repeats_(motif->sequence);
                } else if (tk == "helix3") {
                    result -= params_[1]*num_repeats_(motif->sequence);
                } else if (tk == "helix4") {
                    result -= params_[2]*num_repeats_(motif->sequence);
                } else {
                    result -= params_[3]*num_repeats_(motif->sequence) / ((motif->sequence.size() - 1)/ 2);
                }
            }

            return result;
        }

    private:
        static
        int
        num_repeats_(String const& sequence) {
            assert( sequence.size()%2 == 1 );
            int repeats(0);
            const int helix_size = int((sequence.size() - 1)/2);
            auto curr_pair = String();
            auto next_pair = String();

            for(auto index = 0; index < helix_size - 1; ++index) {
                curr_pair = sequence[index] + *(sequence.rbegin() + index);
                next_pair = sequence[index + 1] + *(sequence.rbegin() + index +1);
                if(curr_pair == next_pair) {
                    repeats += 1;
                }
            }

            return repeats;
        }


    };

}
#endif // RNAMAKE_HELIXREPEATS_H 