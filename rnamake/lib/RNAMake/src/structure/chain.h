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
#include "structure/chain.fwd.h"
#include "structure/residue.h"
#include "structure/residue_type_set.h"

/*
 * Exception for chain
 */
class ChainException : public std::runtime_error {
public:
    /**
     * Standard constructor for ChainException
     * @param   message   Error message for chain
     */
    ChainException(String const & message):
    std::runtime_error(message)
    {}
};

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

class Chain {
public:
    Chain() {}
    Chain(
        ResidueOPs const & residues):
        residues_ ( residues )
    {}
    
    Chain(
        Chain const & c) {
        
        residues_ = ResidueOPs(c.residues_.size());
        int i = 0;
        for (auto const & r : c.residues_) {
            residues_[i] = std::make_shared<Residue>(*r);
            i++;
        }
    }
    
    Chain(
        String const & s,
        ResidueTypeSet const & rts) {
        
        residues_ = ResidueOPs();
        Strings spl = base::split_str_by_delimiter(s, ";");
        for(auto const & r_str : spl) {
            auto r = std::make_shared<Residue>(r_str, rts);
            residues_.push_back(r);
        }
    }
    
    
    ~Chain() {}

public:
    
    inline
    int
    length() const {
        return (int)residues_.size();
    }
    
    inline
    ResidueOP const &
    first() {
        
        if(residues_.size() == 0) {
            throw ChainException("cannot call first there are no "
                                 "residues in chain");
        }
        
        return residues_[0];
    }
    
    inline
    ResidueOP const &
    last() {
        
        if(residues_.size() == 0) {
            throw ChainException("cannot call last there are no "
                                 "residues in chain");
        }
        
        return residues_.back();
    }
    
    inline
    ChainOP
    subchain(int start, int end) {
        if(start < 0) { throw ChainException("start cannot be less then 0"); }
        if(start == end) { throw ChainException("start and end cannot be the same value"); }
        if(end > residues().size()) { throw ChainException("end is greater then chain length"); }
        if(start > end) { throw ChainException("start is greater then end"); }
        
        return ChainOP(new Chain(ResidueOPs(residues_.begin() + start, residues_.begin() + end)));
    }
    
    inline
    ChainOP
    subchain(
        ResidueOP const & r1,
        ResidueOP const & r2) {
        
        if(std::find(residues_.begin(), residues_.end(), r1) == residues_.end()) {
            throw ChainException("starting residue is not part of chain, cannot"
                                 "create subchain");
        
        }
        if(std::find(residues_.begin(), residues_.end(), r2) == residues_.end()) {
            throw ChainException("end residue is not part of chain, cannot"
                                 "create subchain");
        }

        int start = (int)(std::find(residues_.begin(), residues_.end(), r1) - residues_.begin());
        int end   = (int)(std::find(residues_.begin(), residues_.end(), r2) - residues_.begin());
        
        if( start > end) {
            int temp = start;
            start = end;
            end = temp;
        }
        return subchain(start, end);
    }
    
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
    
public: //getters
    
    inline
    ResidueOPs &
    residues() { return residues_; }
    
private:
    ResidueOPs residues_;
};

void
connect_residues_into_chains(
    ResidueOPs & residues,
    ChainOPs & chains);


#endif /* defined(__RNAMake__chain__) */
