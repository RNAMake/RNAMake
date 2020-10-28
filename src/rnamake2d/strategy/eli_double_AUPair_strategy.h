#ifndef RNAMAKE_ELIDOUBLEAUPAIR_H
#define RNAMAKE_ELIDOUBLEAUPAIR_H

#include <regex>
#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class EliDoubleAUPair : public Strategy2D {

    public:
        EliDoubleAUPair() : Strategy2D() {
            params_ = {-3.6974853515624897};
            name_ = "eli_double_AUPair_strategy";
        }

    private:
        Strings
        get_neck_pairs(std::vector<Ints> const& loop_groups, String const& sequence, Ints const& pairmap) {

            auto pairs = Strings{};
            auto pair_indices = Ints{};

            for( auto ii = 0; ii < loop_groups.size(); ++ii ) {
                auto start_i = loop_groups[ii][0];
                auto end_i = loop_groups[ii][(loop_groups[ii]).size() - 1];

                if( start_i - 1 >= 0) {
                    if(pairmap[start_i - 1] < 0) {
                        //print("ERROR something's wrong with RNAELEMENT 1")
                        exit(1);
                    }
                    if(std::find(pair_indices.begin(), pair_indices.end(), start_i - 1) == pair_indices.end()) {
                        auto pair = (
                                String{sequence[std::min(start_i - 1, pairmap[start_i - 1])]}
                                + String{sequence[std::max(start_i - 1, pairmap[start_i - 1])]}
                        );
                        // getting the index and pair
                        std::transform(pair.begin(),pair.end(),pair.end(),toupper);
                        pairs.push_back(
                                std::to_string(std::min(start_i - 1, pairmap[start_i - 1])) + pair
                        );
                        pair_indices.push_back(start_i - 1);
                        pair_indices.push_back(pairmap[start_i - 1]);
                    }
                }
                if( end_i + 1 < sequence.size()) {
                    if (pairmap[end_i + 1] < 0) {
                        //print("ERROR something's wrong with RNAELEMENT 2")
                        exit(1);
                    }
                    if(std::find(pair_indices.begin(), pair_indices.end(), end_i + 1) == pair_indices.end()) {
                        auto pair = (
                                String{sequence[std::min(end_i + 1, pairmap[end_i + 1])]}
                                + String{sequence[std::max(end_i + 1, pairmap[end_i + 1])]}
                        );
                        std::transform(pair.begin(), pair.end(), pair.begin(), toupper);
                        pairs.push_back(std::to_string(std::min(end_i + 1, pairmap[end_i + 1])) + pair);
                        pair_indices.push_back(end_i + 1);
                        pair_indices.push_back(pairmap[end_i + 1]);
                    }
                }
            }

            return pairs;
    }



    public:
        float
        score(Feature2DOP const & feature) override {
            auto neckArea = Strings{};
            auto pairStackGrps = std::vector<Strings>{};
            auto adjPairCounter = 0;
            auto score(100.f);
            for (auto ii = 0; ii < feature->elements.size(); ++ii ) {
                auto elem = feature->elements[ii];
                if( elem.type_ == RNAELEMENT::LOOP) {
                    auto loop_groups = elem.get_loop_groups();
                    auto tmp = get_neck_pairs(loop_groups, feature->sequence, feature->e_pairmap);
                    for(auto jj = 0; jj < tmp.size(); ++jj ){
                        neckArea.push_back(tmp[jj]);
                    }
                }
            }
            for (auto ii = 0; ii < feature->elements.size(); ++ii ) {
                auto elem = feature->elements[ii];
                if( elem.type_ == RNAELEMENT::STACK) {
                    auto lengthOfStack = elem.get_stack_length();
                    auto pairStackRaw = Strings{};
                    auto indexOfPair = Ints{};
                    auto rawIndex = Ints{};

                    for (auto jj = 0; jj < lengthOfStack; ++jj) {
                        rawIndex.push_back(-1);
                    }
                    for (auto jj = 0; jj < lengthOfStack; ++jj) {
                        auto pair = elem.get_pair_from_stack(jj, feature->sequence);

                        if (pair == "UA" or pair == "AU") {
                            indexOfPair.push_back(jj);
                            auto pairContent = (
                                    String{feature->sequence[elem.indices_[jj * 2]]}
                                    + String{feature->sequence[elem.indices_[jj * 2 + 1]]}
                            );
                            std::transform(pairContent.begin(), pairContent.end(), pairContent.begin(), toupper);
                            pairContent = (
                                    std::to_string(std::min(elem.indices_[jj * 2], elem.indices_[jj * 2 + 1]))
                                    + pairContent
                            );
                            pairStackRaw.push_back(pairContent);
                        }
                    }
                    auto counter(0);
                    for (auto jj = 1; jj < indexOfPair.size(); ++jj) {
                        if (indexOfPair[jj] - indexOfPair[jj - 1] == 1) {
                            auto stack = Strings{};
                            stack.push_back(pairStackRaw[jj - 1]);
                            stack.push_back(pairStackRaw[jj]);
                            pairStackGrps.push_back(stack);
                        }
                    }
                }
            }
            auto checker(0);
            const auto digits = std::regex("\\d");
            const auto nondigits = std::regex("\\D");
            //tmp = []
            for(auto kk = 1; kk < pairStackGrps.size(); ++kk )  {
                auto n1 = findall( pairStackGrps[kk - 1][1], digits);
                auto n2 = findall( pairStackGrps[kk][0], digits);
                if(!comp_vectors(n1,n2)) {
                    continue;
                }

                if (
                    comp_vectors(findall( pairStackGrps[kk - 1][0], digits), findall(pairStackGrps[kk - 1][1],digits))
                    and comp_vectors(findall(pairStackGrps[kk][0],digits), findall(pairStackGrps[kk][1],digits))
                    and comp_vectors(findall(pairStackGrps[kk - 1][1],digits),findall(pairStackGrps[kk][0],digits))
                ){
                    //print("this never happen")
                    for(auto ss = 0; ss < neckArea.size(); ++ss)  {
                        if (neckArea[ss] == pairStackGrps[kk - 1][0]
                            or neckArea[ss] == pairStackGrps[kk][1] ) {
                            return 100.f;
                        }
                    }
                }

                pairStackGrps.erase(pairStackGrps.begin() +  kk - 1, pairStackGrps.begin() + kk + 1);
            }
            auto pairStackTuningSame = std::vector<Strings>{};

            for(auto kk = 0; kk < pairStackGrps.size(); ++kk) {
                for(auto mm  = 1; mm < pairStackGrps[kk].size(); ++mm ) {

                    auto a = findall(pairStackGrps[kk][mm - 1], nondigits);
                    auto b = findall(pairStackGrps[kk][mm], nondigits);
                    if ( comp_vectors(a,b)) {
                        auto stack = Strings{};
                        stack.push_back(pairStackGrps[kk][mm - 1]);
                        stack.push_back(pairStackGrps[kk][mm]);
                        pairStackTuningSame.push_back(stack);
                    }
                }
            }
            adjPairCounter = pairStackTuningSame.size();
            for(auto kk = 0; kk < pairStackTuningSame.size(); ++kk) {
                auto numOfNeck(0);
                for( auto mm  = 0; mm <  pairStackTuningSame[kk].size(); ++mm ) {
                    for(auto yy  = 0; yy < neckArea.size(); ++yy) {
                        if ( pairStackTuningSame[kk][mm] == neckArea[yy] ) {
                            ++numOfNeck;
                        }
                    }
                    }

                if( numOfNeck > 0) {
                    adjPairCounter = adjPairCounter - 1;
                }
            }
            if (adjPairCounter > 0 ) {
                score += (adjPairCounter - 1) * params_[0];
            }
            return score;
            }
    };

}
#endif // RNAMAKE_ELIDOUBLEAUPAIR_H 