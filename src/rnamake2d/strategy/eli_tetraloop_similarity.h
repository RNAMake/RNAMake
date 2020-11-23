#ifndef RNAMAKE_ELITETRALOOPSIMILARITY_H
#define RNAMAKE_ELITETRALOOPSIMILARITY_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class EliTetraloopSimilarity : public Strategy2D {
    private:

    public:
        EliTetraloopSimilarity() : Strategy2D() {
            params_ = {4.f, 4.f, 0.3f};
            name_ = "eli_tetraloop_similarity";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            const auto& sequence = feature->sequence;
            const auto& elements = feature->elements;
            const auto& pairmap = feature->pairmap;

            auto penalty(0);
            auto count(0);
            auto tetraloop_idx = Ints();
            auto energies = std::vector<float>();
            for(int ii = 0; ii < elements.size(); ++ii) {
                if(elements[ii].type_ == RNAELEMENT::LOOP
                    && elements[ii].children_.size() == 0
                    && elements[ii].parent_
                    && elements[ii].indices_.size() == 4
                ) {
                    count += 1;
                    tetraloop_idx.push_back(ii);
                    energies.push_back(elements[ii].score_);
                }
            }
            auto score(0.f);
            if(count == 4) {
                 std::sort(energies.begin(), energies.end());
                if (energies[3] - energies[0] < params_[2]) {
                    score += 5.f;
                } else if (
                        energies[3] - energies[2] < params_[2]
                        and energies[1] - energies[0] < params_[2]
                ) {
                    score += 5.f;
                } else if (
                        energies[3] - energies[1] < params_[2]
                        or energies[2] - energies[0] < params_[2]
                ) {
                    score += 3.f;
                } else if (
                        energies[3] - energies[2] < params_[2]
                        or energies[2] - energies[1] < params_[2]
                        or energies[1] - energies[0] < params_[2]
                ) {
                    score += 2.f;
                }
                return score*20.f;
            } else if (count == 3) {
                std::sort(energies.begin(), energies.end());
                if( energies[2] - energies[0] < params_[2]) {
                    score += 5.f;
                } else if (
                        energies[2] - energies[1] < params_[2]
                        or energies[1] - energies[0] < params_[2]
                ) {
                    score += 3.f;
                }
                return score * 20.f;
            } else {
                return UNSCORABLE;
            }
        }
    };

}
#endif // RNAMAKE_ELITETRALOOPSIMILARITY_H 