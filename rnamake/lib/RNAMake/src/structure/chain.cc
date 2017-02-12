//
//  chain.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <map>
#include "chain.h"

ChainOP
Chain::subchain(int start, int end) {
    if (start < 0)              { throw ChainException("start cannot be less then 0"); }
    if (start == end)           { throw ChainException("start and end cannot be the same value"); }
    if (end > residues_.size()) { throw ChainException("end is greater then chain length"); }
    if (start > end)            { throw ChainException("start is greater then end"); }

    return std::make_shared<Chain>(ResidueOPs(residues_.begin() + start, residues_.begin() + end));
}

ChainOP
Chain::subchain(
        ResidueOP const & r1,
        ResidueOP const & r2) {

    if (std::find(residues_.begin(), residues_.end(), r1) == residues_.end()) {
        throw ChainException("starting residue is not part of chain, cannot"
                                     "create subchain");

    }
    if (std::find(residues_.begin(), residues_.end(), r2) == residues_.end()) {
        throw ChainException("end residue is not part of chain, cannot"
                                     "create subchain");
    }

    int start = (int) (std::find(residues_.begin(), residues_.end(), r1) - residues_.begin());
    int end = (int) (std::find(residues_.begin(), residues_.end(), r2) - residues_.begin());

    if (start > end) {
        int temp = start;
        start = end;
        end = temp;
    }
    return subchain(start, end);
}

String
Chain::to_str() const {
    String s;
    for (auto const & r : residues_) {
        s += r->to_str() + ";";
    }
    return s;
}

String
Chain::to_pdb_str(
        int & acount,
        int rnum,
        String const & chain_id) const {

    if (rnum == -1) {
        rnum = residues_[0]->num();
    }

    String s;
    for (auto const & r : residues_) {
        s += r->to_pdb_str(acount, rnum, chain_id);
        rnum += 1;
    }
    return s;
}

void
Chain::to_pdb(
        String const fname,
        int rnum,
        String const & chain_id) const {
    std::ofstream out;
    out.open(fname.c_str());
    int i = 1;
    String s = to_pdb_str(i, rnum, chain_id);
    out << s << std::endl;
    out.close();
}

void
connect_residues_into_chains(
        ResidueOPs & residues,
        ChainOPs & chains) {

    ResidueOP current;
    ResidueOPs current_chain_res;
    int five_prime_end = 1;
    int found = 1;
    while (true) {
        for (auto & r1 : residues) {
            five_prime_end = 1;
            for (auto & r2 : residues) {
                if (r1->connected_to(*r2) == -1) {
                    five_prime_end = 0;
                    break;
                }
            }
            if (five_prime_end) {
                current = r1;
                break;
            }
        }
        if (!five_prime_end) { break; }
        residues.erase(std::remove(residues.begin(), residues.end(), current), residues.end());
        current_chain_res = ResidueOPs();
        found = 1;
        while (found) {
            current_chain_res.push_back(current);
            found = 0;
            for (auto & r : residues) {
                if (current->connected_to(*r) == 1) {
                    current = r;
                    found = 1;
                    break;
                }
            }
            if (found) {
                residues.erase(std::remove(residues.begin(), residues.end(), current), residues.end());
            } else {
                chains.push_back(ChainOP(new Chain(current_chain_res)));
            }
        }
        if (residues.size() == 0) {
            break;
        }
    }

}

ResiduesandChainCuts
get_chain_cuts(ResidueOPs const & residues) {

    ResidueOP current;
    ResidueOPs current_chain_res;
    ResidueOPs all_res;
    Ints chain_cuts;
    auto seen = std::map<ResidueOP, int>();
    int five_prime_end = 1;
    int found = 1;
    while (true) {
        for (auto & r1 : residues) {
            if (seen.find(r1) != seen.end()) { continue; }
            five_prime_end = 1;
            for (auto & r2 : residues) {
                if (seen.find(r2) != seen.end()) { continue; }
                if (r1->connected_to(*r2) == -1) {
                    five_prime_end = 0;
                    break;
                }
            }
            if (five_prime_end) {
                current = r1;
                break;
            }
        }
        if (!five_prime_end) { break; }
        seen[current] = 1;
        current_chain_res = ResidueOPs();
        found = 1;
        while (found) {
            current_chain_res.push_back(current);
            found = 0;
            for (auto & r : residues) {
                if (seen.find(r) != seen.end()) { continue; }
                if (current->connected_to(*r) == 1) {
                    current = r;
                    found = 1;
                    break;
                }
            }
            if (found) {
                seen[current] = 1;
            } else {
                for(auto const & r : current_chain_res) {
                    all_res.push_back(r);
                }
                chain_cuts.push_back((int)all_res.size());

                break;
            }
        }
        if(seen.size() == residues.size()) { break; }
    }

    chain_cuts.push_back((int)all_res.size());

    return ResiduesandChainCuts { all_res, chain_cuts};

}

bool
are_chains_equal(
        ChainOP const & c1,
        ChainOP const & c2,
        int check_uuids) {

    if (c1->length() != c2->length()) { return false; }

    auto result = 0;
    for (int i = 0; i < c1->length(); i++) {
        result = are_residues_equal(c1->get_residue(i), c2->get_residue(i), check_uuids);
        if (!result) { return false; }
    }
    return true;
}
