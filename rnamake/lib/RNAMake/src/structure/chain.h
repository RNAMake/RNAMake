//
//  chain.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__chain__
#define __RNAMake__chain__

#include <stdio.h>
#include <algorithm>

//RNAMake Headers
#include "base/types.h"
#include "primitives/chain.h"
#include "motif_state/chain.h"
#include "structure/chain.fwd.h"
#include "structure/residue.h"
#include "structure/residue_type_set.h"

/**
 * Stored chain information from pdb file. Stores all residues in chain.
 * Implementation is designed to be extremely lightweight. To connect residues
 * into chains it highly adviced that you use function: connect_residues_into_chains
 *
 * @code
 *  //grabbing example instance
 *  #include "instances/structure_instances.hpp" 
 *  auto c = instances::chain();
 *
 *  std::cout << c->first() << std::endl;
 *  //OUTPUT: <Residue('G103 chain A')>
 *
 *  std::cout << c->last() << std::endl;
 *  //OUTPUT: <Residue('C260 chain A')>
 *
 *  std::cout << c->length() << std::endl;
 *  //OUTPUT: 157
 *
 *  auto cs = c->subchain(1, 10)
 *  std::cout << cs.length() << std::endl;
 *  //OUTPUT: 9
 *
 *  std::cout << cs.first() << std::endl;
 *  //OUTPUT: <Residue('A104 chain A')>
 *
 *  auto cs2 = c->subchain(start_res=c->residues()[10], end_res=c->residues()[15])
 *  std::cout << cs2.lengt() << std::endl;
 *  //OUTPUT: 6
 *
 *  c->to_pdb_str()
 *  //OUTPUT:
 *  ATOM      1 O5'  G   A 103     -26.469 -47.756  84.669  1.00  0.00
 *  ATOM      2 C5'  G   A 103     -25.050 -47.579  84.775  1.00  0.00
 *  ATOM      3 C4'  G   A 103     -24.521 -48.156  86.068  1.00  0.00
 *  ATOM      4 O4'  G   A 103     -24.861 -49.568  86.118  1.00  0.00
 *  ATOM      5 C3'  G   A 103     -23.009 -48.119  86.281  1.00  0.00
 *  ATOM      6 O3'  G   A 103     -22.548 -46.872  86.808  1.00  0.00
 *  ATOM      7 C1'  G   A 103     -23.806 -50.289  86.732  1.00  0.00
 *  ATOM      8 C2'  G   A 103     -22.812 -49.259  87.269  1.00  0.00
 *  .
 *  .
 *  .
 */

class Chain : public primitives::Chain<Residue> {
public:
    Chain(
            ResidueOPs const & residues) :
            primitives::Chain<Residue>(residues) {}

    Chain(
            Chain const & c,
            int new_uuid = 0) :
            primitives::Chain<Residue>() {

        residues_ = ResidueOPs(c.residues_.size());
        int i = 0;
        for (auto const & r : c.residues_) {
            residues_[i] = std::make_shared<Residue>(*r, new_uuid);
            i++;
        }
    }

    Chain(
            String const & s,
            ResidueTypeSet const & rts) :
            primitives::Chain<Residue>() {

        residues_ = ResidueOPs();
        Strings spl = split_str_by_delimiter(s, ";");
        for (auto const & r_str : spl) {
            auto r = std::make_shared<Residue>(r_str, rts);
            residues_.push_back(r);
        }
    }

    ~Chain() {}

public:

    ChainOP
    subchain(int, int);

    ChainOP
    subchain(
            ResidueOP const &,
            ResidueOP const &);

    String
    to_str() const;


    String
    to_pdb_str(
            int &,
            int,
            String const &) const;

    inline
    String
    to_pdb_str(int & acount) const {
        return to_pdb_str(acount, -1, "");
    }

    void
    to_pdb(
            String const,
            int,
            String const &) const;

    inline
    void
    to_pdb(String const & fname) {
        return to_pdb(fname, -1, "");
    }

    state::ChainOP
    get_state() {
        auto residues = state::ResidueOPs();
        for(auto const & r : residues_) {
            residues.push_back(r->get_state());
        }
        return std::make_shared<state::Chain>(residues);
    }

public:

    inline
    void
    move(Point const & p) {
        for(auto & r : residues_) { r->move(p); }
    }

    inline
    void
    transform(Transform const & t) {
        for(auto & r : residues_) { r->transform(t); }
    }

    inline
    void
    fast_transform(
            Matrix const & r,
            Vector const & t) {
        for(auto & res : residues_) { res->fast_transform(r, t); }
    }


};

struct ResiduesandChainCuts {
    ResidueOPs residues;
    Ints chain_cuts;
};

void
connect_residues_into_chains(
        ResidueOPs &,
        ChainOPs &);

ResiduesandChainCuts
get_chain_cuts(ResidueOPs const &);

bool
are_chains_equal(
        ChainOP const & c1,
        ChainOP const & c2,
        int check_uuids = 1);

#endif /* defined(__RNAMake__chain__) */
