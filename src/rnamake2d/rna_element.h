#ifndef RNAMAKE_RNA_ELEMENT_H
#define RNAMAKE_RNA_ELEMENT_H

#include <regex>
#include <algorithm>
#include <vector>

#include <base/types.h>
#include <plog/Log.h>

namespace rnamake2d {
    constexpr auto UNSCORABLE = -99999;
    constexpr auto SHIFT_LIMIT = 3;
    constexpr auto PAIR_TYPES = 3;
    constexpr auto TOLERANCE = 1e-3f;

    enum RNAELEMENT {
        LOOP,
        STACK
    };
    struct RNAElement {
        Ints indices_ = {};
        int type_ = -1;
        int branching_stacks_ = 0;
        float score_ = 0.f;
        RNAElement* parent_ = nullptr;
        std::vector<RNAElement*> children_ = {};
        Ints quad_scores_ = {};

        std::vector<Ints>
        get_loop_groups() const {
            auto groups = std::vector<Ints>{} ;
            if (type_ != RNAELEMENT::LOOP) {
                return groups;
            }
            auto last_index = -999;
            auto last_group = Ints{};
            for(auto ii = 0; ii < indices_.size(); ++ii) {
                if (indices_[ii] != last_index + 1 and last_index >= 0) {
                    groups.push_back(last_group);
                    last_group.clear();
                }
                last_group.push_back(indices_[ii]);
                last_index = indices_[ii];
            }

            if(!last_group.empty()) {
                groups.push_back(last_group);
            }

            return groups;
        }
        int
        get_stack_length() const  {
            if (type_ != RNAELEMENT::STACK) {
                return 0;
            }
            return int(indices_.size()/ 2);
        }

        String
        get_pair_from_stack(int  pair_index,String const& sequence) const {
            auto pair = String{sequence[indices_[int(pair_index) * 2]]};
            pair += sequence[indices_[int(pair_index) * 2 + 1]];
            std::transform(pair.begin(), pair.end(), pair.begin(), toupper);
            return  pair;
        }

        Strings
        get_loop_closing_pairs(String const&  sequence, Ints const & pairmap) const {
            auto pairs =  Strings{};
            if( type_ != RNAELEMENT::LOOP) {
                return pairs;
            }
            auto pair_indices =  std::vector<Ints>{};
            if (parent_ != nullptr) {
                const auto& parent = parent_;
                if (parent->type_ != RNAELEMENT::STACK) {
                    //print("ERROR Loop's parent is not a stack")
                    exit(1);
                }
                const auto npi = parent->indices_.size();
                String pair = String{sequence[parent->indices_[npi - 2]]} + String{sequence[parent->indices_[npi - 1]]};
                pairs.push_back(pair);
            }
            for( auto ii = 0; ii < children_.size(); ++ii ) {
                const auto& child = children_[ii];
                if(child->type_ != RNAELEMENT::STACK) {
                    //print("ERROR Loop's child is not a stack")
                    exit(0);
                }
                const auto nci = child->indices_.size();
                pairs.emplace_back(String{sequence[child->indices_[0]]} + String{sequence[child->indices_[1]]});
            }
            return pairs;
        }
    }; // class RNAElement

    using RNAElems = std::vector<RNAElement>;

    Ints
    get_pairmap_from_secstruct(String const& );

    RNAElems
    get_rna_elemnts_from_secstruct(String const& );

    void
    get_rna_elements_from_secstruct_recursive(Ints const&, int, int, RNAElems&, int, int, RNAElement*);


    Strings
    findall(String const& , std::regex const&) ;

    template<typename T>
    bool
    comp_vectors(std::vector<T> const& v1, std::vector<T> const& v2 ) {
        if(v1.size() != v2.size())  {
            return false;
        } else {
            for(auto ii = 0; ii < v1.size(); ++ii) {
                if(v1[ii] != v2[ii]) {
                    return false;
                }
            }
            return true;
        }
    }
} // namespace rnamake2d

#endif // RNAMAKE_RNA_ELEMENT_H