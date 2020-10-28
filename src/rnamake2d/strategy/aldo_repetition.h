#ifndef RNAMAKE_ALDOREPETITION_H
#define RNAMAKE_ALDOREPETITION_H

#include <unordered_set>
#include <unordered_map>

#include <base/types.h>
#include <base/exception.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class AldoRepetition : public Strategy2D {
    private:

    public:
        AldoRepetition() : Strategy2D() {
            params_ = { 9.30186684,
                        1.05847095,
                        14.18676777,
                        -25.33226257};
            name_ = "AldoRepetition";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto repeat_2(0), repeat_3(0), repeat_4(0), repeat_5(0);

            auto sub_seq_dict = std::unordered_map<String,int>{};

            const auto length = feature->structure.size();

            for(auto start_index = 0; start_index < length; ++start_index) {
                for(auto seq_len = 2; seq_len < 6; ++seq_len) {
                    if(start_index + seq_len >= length) {
                        break;
                    }
                    const auto sub_seq = feature->sequence.substr(start_index, seq_len);
                    sub_seq_dict[sub_seq]++;
                }
            }

            for(const auto& pr : sub_seq_dict) {
                const auto seq_len = pr.first.size();
                if(seq_len < 2) {
                    throw base::RNAMakeException("Repeat sequence length should not be less than 2");
                } else if (seq_len == 2 ) {
                    ++repeat_2;
                } else if (seq_len == 3) {
                    ++repeat_3;
                } else if (seq_len == 4) {
                    ++repeat_4;
                } else if (seq_len == 5) {
                    ++repeat_5;
                }
            }

            return 100.f - repeat_2*params_[0] - repeat_3*params_[1] - repeat_4*params_[2] - repeat_5*params_[3];
        }
    };

}
#endif // RNAMAKE_ALDOREPETITION_H 