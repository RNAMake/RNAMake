#ifndef RNAMAKE_BAD_HELIX_H
#define RNAMAKE_BAD_HELIX_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class BadHelix4 : public Strategy2D  {
    private:
        std::unordered_set<String> bad_ = {
                "ACUG&UAGU",
                "CAGA&UCUG",
                "CUAU&AUCC",
                "GAGA&UCUC",
                "GCUU&GAGU",
                "GGGA&UCCU",
                "UUAU&AUGA",
                "UUCU&AGAA",
        };

    public:
        BadHelix4() : Strategy2D() {
                params_ = {12.20549451};
                name_ = "BadHelix4";
        }

    public:
        float
        score(Feature2DOP const & feature ) override {
            auto result(100.f);

            for(const auto& motif : feature->motifs) {
                if(motif->mtype() == util::MotifType::HELIX
                    && bad_.find(motif->sequence) != bad_.end()
                )  {
                    result -= params_[1];
                }
            }

            return result;
        }
    };
}

#endif // RNAMAKE_BAD_HELIX_H