#ifndef RNAMAKE_ELITETRALOOPBLUES_H
#define RNAMAKE_ELITETRALOOPBLUES_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class EliTetraloopBlues : public Strategy2D {
    private:

    public:
        EliTetraloopBlues() : Strategy2D() {
            params_ = {4.f};
            name_ = "eli_tetraloop_blues";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            const auto& sequence = feature->sequence;
            const auto& elements = feature->elements;
            const auto& pairmap = feature->e_pairmap;

            auto penalty(0), count(0);

            for(auto ii = 0; ii < elements.size(); ++ii) {
                if(elements[ii].type_ == RNAELEMENT::LOOP
                    and elements[ii].children_.size() == 0
                    and elements[ii].parent_
                    and elements[ii].indices_.size() == 4 ) {
                        count++;
                        const auto indices = elements[ii].indices_;
                        auto u_count(0);
                        for(auto jj = 0; jj < indices.size(); ++jj) {
                            if(sequence[indices[jj]] == 'U') {
                                u_count++;
                            }
                        }

                        const auto& parent = elements[ii].parent_;
                        const auto npi = parent->indices_.size();
                        for(auto kk = 1; kk < 3; ++kk) {
                            if(sequence[parent->indices_[npi -kk]] == 'U') {
                                u_count++;
                            }
                        }

                        if(u_count > 2 ) {
                            penalty += (u_count - 2);
                        }
                    }
                }

            if(count == 0)  {
                return UNSCORABLE;
            }
            return 100.f - float(penalty)*params_[0];
        }
    };

}
#endif // RNAMAKE_ELITETRALOOPBLUES_H 