#ifndef RNAMAKE_ELILEGALPLACEMENTOFGUPAIRS_H
#define RNAMAKE_ELILEGALPLACEMENTOFGUPAIRS_H

#include <unordered_set>

#include <base/types.h>
#include <rnamake2d/strategy2d.h>

namespace rnamake2d {
    class EliLegalPlacementOfGUPairs : public Strategy2D {
    private:

    public:
        EliLegalPlacementOfGUPairs() : Strategy2D() {
            params_ = {-0.5600628139854402,
                       1.0642537595877468,
                       0.21755765796198936,
                       1.6772074364078484,
                       2.473752807629996,
                       2.0126349545807547,
                       -3.828471068768436,
                       -1.5418531268055973,
                       1.0938302617479152,
                       -0.9704340785521139,
                       -2.1808586124069977,
                       -2.5046531589390852,
                       -1.680199359651886,
                       1.7905500957241554
                        };
            name_ = "eli_legal_placement_of_GUpairs";
        }

    private:
        static
        std::vector<Strings>
        get_original_pairs(RNAElement const & elem, String const& sequence) {
            auto pairs = Strings{};
            auto pairs_sep = std::vector<Strings>{};
            for(auto i = 0; i < elem.get_stack_length(); ++i) {
                const auto pair = elem.get_pair_from_stack(i, sequence);
                pairs.push_back(pair);
            }
            auto tmp(0);
            auto p = Strings{};
            for(auto i = 0; i < pairs.size(); ++i) {
                const auto& pair = pairs[i];
                if(pair == "GU" || pair == "UG") {
                    for(auto j = tmp; j < i; ++j)  {
                        p.push_back(pairs[i]);
                    }

                    if(p.size() > 0) {
                        pairs_sep.push_back(p);
                    } else if (tmp == 0) {
                        pairs_sep.push_back(p);
                    }
                    tmp = i + 1;
                }
            }
            if(!p.empty())  {
                p.clear();
            }

            if(tmp != 0) {
                for(auto j = tmp; j < pairs.size(); ++j) {
                    p.push_back(pairs[j]);
                }
                pairs_sep.push_back(p);
            }

            return pairs_sep;
        }

    private:
        static
        std::vector<Strings>
        get_filted_pairs(RNAElement const& elem, String const & sequence ) {
            auto pairs_sep = get_original_pairs(elem, sequence);
            auto pairs_sep_filted = std::vector<Strings>{};

            for(auto ii = 0; ii < pairs_sep.size(); ++ii) {
                auto tmp = Strings{};
                for(auto jj = 0; jj < pairs_sep[ii].size(); ++jj) {
                    if(pairs_sep[ii][jj] == "GC" || pairs_sep[ii][jj] == "CG") {
                        tmp.push_back(pairs_sep[ii][jj]);
                    }
                }
                pairs_sep_filted.push_back(tmp);
            }
            return pairs_sep_filted;
        }

    private:
        static
        int
        get_num_twistGC(Strings const& groupA, Strings const& groupB) {
            auto result(0);
            for(const auto& pairA : groupA) {
                for(const auto& pairB : groupB) {
                    if(pairA != pairB) {
                        ++result;
                    }
                }
            }
            return result;
        }

    private:
        static
        int
        get_num_sameGC(Strings const& groupA, Strings const& groupB) {
            auto result(0);
            for(const auto& pairA : groupA) {
                for(const auto& pairB : groupB) {
                    if(pairA == pairB) {
                        ++result;
                    }
                }
            }
            return result;
        }

    private:
        static
        Ints
        two_GC(RNAElement const & elem, String const& sequence ) {

            auto numOfTwist(0);
            auto numOfSame(0);
            auto result = Ints{};
            auto pairs_sep_filted = get_filted_pairs(elem, sequence);

            if(pairs_sep_filted.empty()) {
                return result;
            }
            for(int ii = 0; ii < pairs_sep_filted.size() - 1; ++ii) {
                const auto pair = pairs_sep_filted[ii];
                for(auto jj = ii + 1; jj < pairs_sep_filted.size(); ++jj) {
                    const auto compare_pair = pairs_sep_filted[jj];

                    numOfTwist += get_num_twistGC(pair, compare_pair);
                    numOfSame += get_num_sameGC(pair, compare_pair);
                }
            }

            result.push_back(numOfTwist);
            result.push_back(numOfSame);

            return result;
        }

    private:
        static
        int
        besideGC(RNAElement const & elem, String const& sequence) {
            auto result(0);
            const auto pairs = get_filted_pairs(elem, sequence);

            if(pairs.size() < 2) {
                return result;
            }

            const auto tmp = pairs[0].size();

            if(tmp != 0 && (pairs[0][tmp - 1] == "GC" && pairs[0][tmp-1] == "CG" )) {
                ++result;
            }

            if(pairs.size() > 0) {
                const auto tmp = pairs.size() - 1;
                if(pairs[tmp].size() > 0 && (pairs[tmp][0] == "GC" || pairs[tmp][0] == "CG"))  {
                    ++result;
                }
            }

            if(pairs.size() > 2) {
                for(auto ii = 0; ii < pairs.size() - 1; ++ii) {
                    const auto tmp = pairs[ii].size();

                    if(tmp > 0 && (pairs[ii][tmp -1] == "GC" || pairs[ii][tmp - 1] == "CG")) {
                        ++result;
                    }

                    if(pairs[ii].size() > 0 &&  (pairs[ii][0] == "GC" || pairs[ii][0] == "CG")) {
                        ++result;
                    }
                }
            }

            return result;
        }

    private:
        static
        int
        get_besidingGU_pairs(RNAElement const & elem, String const & sequence) {
            auto result(0);
            auto pairs = Strings{};

            for(auto i = 0; i < elem.get_stack_length(); ++i) {
                const auto pair = elem.get_pair_from_stack(i, sequence);
                pairs.push_back(pair);
            }

            for(auto i = 0; i < pairs.size() - 1 ; ++i) {
                if((pairs[i] == "GU" || pairs[i] == "UG") &&
                    (pairs[i+1] == "GU" || pairs[i+1] == "UG")) {
                    ++result;
                }
            }
            return result;
        }

    public:
        float
        score(Feature2DOP const & feature) override {
            auto result(80.f);
            auto noGU(true);
            for(auto ii = 0; ii < feature->elements.size(); ++ii) {
                const auto& elem = feature->elements[ii];
                if(elem.type_ == RNAELEMENT::STACK ){
                    if(get_original_pairs(elem, feature->sequence).size() != 0)  {
                        noGU = false;
                        break;
                    }
                }
            }


            if(noGU) {
                result += params_[0];
            }

            for(auto ii = 0; ii < feature->elements.size(); ++ii)  {
                const auto& elem = feature->elements[ii];
                if(elem.type_ == RNAELEMENT::STACK) {
                    const auto numBetween = two_GC(elem, feature->sequence);
                    if(numBetween.empty()) {
                        continue;
                    }
                    result += numBetween[0]*params_[1];
                    result += numBetween[1]*params_[2];
                    const auto numBeside = besideGC(elem, feature->sequence);
                    if(numBeside != 0) {
                        result += numBeside*params_[3];
                    } else {
                        result += params_[4];
                    }
                }
            }


            auto besidingGU(0);

            for(auto ii = 0; ii < feature->elements.size(); ++ii) {
                const auto elem = feature->elements[ii];
                if(elem.type_ == RNAELEMENT::STACK) {
                    besidingGU += get_besidingGU_pairs(elem, feature->sequence);
                }
            }

            if(besidingGU == 1) {
                result += params_[5];
            } else if (besidingGU > 1) {
                result += float(besidingGU - 1) * params_[6];
            }

            for(auto ii = 0; ii < feature->elements.size(); ++ii) {
                const auto elem = feature->elements[ii];

                if(elem.type_ == RNAELEMENT::STACK) {
                    const auto pairs = get_original_pairs(elem, feature->sequence);
                    if(pairs.size() == 0) {
                        continue;
                    }
                    if(pairs[0].size() == 0 or pairs[pairs.size() - 1].size() == 0) {
                        result += params_[7];
                    }
                }
            }

            auto indeOfNeck = -1;
            auto first_nuc_with_pair = -1;
            auto last_nuc_with_pair = -1;

            for(auto i = 0; i < feature->sequence.size(); ++i) {
                if(i != -1) {
                    first_nuc_with_pair = i;
                    last_nuc_with_pair = feature->pairmap[i];
                }
            }

            RNAElement neckArea;
            for(auto ii = 0; ii < feature->elements.size(); ++ii) {
                const auto elem = feature->elements[ii];
                if(elem.type_ == RNAELEMENT::STACK) {
                    for(auto jj = 0 ; jj < elem.get_stack_length(); ++jj) {
                        if(elem.indices_[jj] == first_nuc_with_pair) {
                            break;
                        }
                    }
                    neckArea = elem;
                    break;
                }
            }

            const auto numBetween = two_GC(neckArea, feature->sequence);
            if(!numBetween.empty())  {
                result += numBetween[0] * params_[8];
                result += numBetween[1] * params_[9];
            }

            const auto numBeside = besideGC(neckArea, feature->sequence);
            if(numBeside != 0) {
                result += numBeside*params_[10];
            } else {
                result += params_[11];
            }

            auto pairs = Strings{};
            auto tmp(0);
            auto isSameTuning(true);

            for(auto i = 0; i < neckArea.get_stack_length(); ++i) {
                const auto pair = neckArea.get_pair_from_stack(i, feature->sequence);
                pairs.push_back(pair);
            }
            if(!pairs.empty()) {
                for (int i = 0; i < pairs.size() - 1; ++i) {
                    if ((pairs[i] == "GU" || pairs[i] == "UG")
                            and (pairs[i + 1] == "GU" || pairs[i + 1] == "UG") ) {
                        ++tmp;
                        if (pairs[i] != pairs[i + 1]) {
                            isSameTuning = false;
                        }
                    }
                }
            }

            if(tmp > 1 || isSameTuning) {
                result += params_[12];
            }

            if(!pairs.empty()) {
                if (pairs[0] == "GU"
                    || pairs[0] == "UG"
                    || *pairs.rbegin() == "GU"
                    || *pairs.rbegin() == "UG"
                        ) {
                    result += params_[13];
                }
            }

            return result;
        }
    };

}
#endif // RNAMAKE_ELILEGALPLACEMENTOFGUPAIRS_H 