#ifndef RNAMAKE_HELIX_2_H
#define RNAMAKE_HELIX_2_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class Helix2 : public Strategy2D  {
    private:
    std::unordered_set<String> good_ = {
            "AC&GU", "AG&CU", "CC&GG", "CG&CG", "GC&GC", "GG&CC", "GU&AC"};
    public:
        Helix2() : Strategy2D() {
            params_ = {9.20041156e+01, -7.82690044e-03, -6.88522970e-01,  2.68064998e+00};
            name_ = "Helix2";
        }

    public:
        float
        score(Feature2DOP const & feature ) override {
            auto result(params_[0]);

            for(const auto& motif : feature->motifs)  {
                if(motif->mtype() != util::MotifType::HELIX ||
                motif->buffer_ != 2) {
                    continue;
                }

                if(good_.find(motif->sequence) != good_.end()) {
                    if(motif->sequence.find('A') != String::npos) {
                       result += params_[2];
                    } else {
                        result += params_[1];
                    }
                } else {
                    result -= params_[3];
                }
            }
            return result;
        }
    };
}

#endif// RNAMAKE_HELIX_2_H
