#ifndef RNAMAKE_SINGLESTRAND2_H
#define RNAMAKE_SINGLESTRAND2_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class SingleStrand2 : public Strategy2D {
    private:
        std::unordered_set<String> good_ = {
                "AA"
        };
    public:
        SingleStrand2() : Strategy2D() {
            params_ = {1.};
            name_ = "SingleStrand2";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);
            for(const auto& ss : feature->motifs) {
                if(ss->mtype() != util::MotifType::SSTRAND || ss->sequence.size() != 3) {
                    continue;
                }
                if(good_.find(ss->sequence) != good_.end()) {
                    result -= params_[0];
                }
            }
            return result;
        }
    };

}
#endif // RNAMAKE_SINGLESTRAND2_H 