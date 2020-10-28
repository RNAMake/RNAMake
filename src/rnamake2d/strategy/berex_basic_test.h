#ifndef RNAMAKE_BEREXBASICTEST_H
#define RNAMAKE_BEREXBASICTEST_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class BerexBasicTest : public Strategy2D {
    private:

    public:
        BerexBasicTest() : Strategy2D() {
            params_ = {0.1949575184053913,
                       46.2935845771954,
                       0.1161491442538988,
                       36.58439589059461,
                       0.18909543838424497,
                       48.543169273442246,
                       -219.59256786408753,
                       -38.07400199637601,
                       -0.0015152832056013763,
                       69.64452390684735,
                       150.15689003443046,
                       1.9660209419220278};
            name_ = "berex_basic_test";
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(100.f);
            auto a_count(0), c_count(0), g_count(0), u_count(0);

            const auto& seq = feature->sequence;
            const auto seq_len = seq.size();

            for(auto ii = 0; ii < seq_len; ++ii) {
                if(seq[ii] == 'G') {
                    ++g_count;
                } else if (seq[ii] == 'C') {
                    ++c_count;
                } else if (seq[ii] == 'A') {
                    ++a_count;
                } else {
                    ++u_count;
                }
            }

            result -= std::abs(float(g_count) / float(seq_len) - params_[0])*params_[1];
            result -= std::abs(float(u_count) / float(seq_len) - params_[2])*params_[3];
            result -= std::abs(float(c_count) / float(seq_len) - params_[4])*params_[5];
            const auto fe = feature->fe;

            if(fe < params_[6]) {
                result -= std::abs(fe - params_[6]) *params_[8];
            } else if (fe > params_[7]) {
                result -= std::abs(fe - params_[7]) *params_[8];
            }

            const auto mp = feature->meltpoint;
            if(mp < params_[9]) {
                result -= std::abs(mp - params_[9]) * params_[11];
            } else if (mp > params_[10]) {
                result -= std::abs(mp - params_[10]) * params_[11];
            }

            return result;
        }
    };

}
#endif // RNAMAKE_BEREXBASICTEST_H 