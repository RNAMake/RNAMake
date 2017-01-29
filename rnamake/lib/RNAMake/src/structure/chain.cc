//
//  chain.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

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
