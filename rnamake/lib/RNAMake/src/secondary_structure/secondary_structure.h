//
//  secondary_structure.h
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__secondary_structure__
#define __RNAMake__secondary_structure__

#include <stdio.h>
#include <memory.h>
#include <map>

#include "base/types.h"
#include "util/uuid.h"

#include "secondary_structure/chain.h"
#include "secondary_structure/motif.h"

namespace sstruct {

    
class SecondaryStructure : public Motif {
public:
    
    SecondaryStructure() {}
    
    SecondaryStructure(
        String const &,
        String const &);
    
    SecondaryStructure(
        ChainOPs const &);
    
    ~SecondaryStructure() {}
    
    
public:
    
    inline
    SecondaryStructure
    copy() {
        ChainOPs new_chains;
        for(auto const & c : chains_) {
            new_chains.push_back(std::make_shared<Chain>(c->copy()));
        }
        
        BasepairOPs bps, ends;
        auto n_ss = SecondaryStructure(new_chains);
        for(auto const & bp : basepairs_) {
            auto new_bp = Basepair(n_ss.get_residue(bp->res1()->uuid()),
                                   n_ss.get_residue(bp->res2()->uuid()),
                                   bp->uuid());
            auto new_bp_pointer = std::make_shared<Basepair>(new_bp);
            bps.push_back(new_bp_pointer);
        }
        for(auto const & end : ends_) {
            int i = (int)(std::find(basepairs_.begin(), basepairs_.end(), end) - basepairs_.begin());
            ends.push_back(bps[i]);
        }
        n_ss.basepairs(bps);
        n_ss.ends(ends);
        
        std::map<String, MotifOPs> motifs;
        motifs["ALL"] = MotifOPs();
        for(auto const & m : motifs_["ALL"]) {
            ChainOPs m_new_chains;
            for(auto const & c : m->chains()) {
                ResidueOPs new_res;
                for(auto const & r : c->residues()) {
                    auto n_r = n_ss.get_residue(r->uuid());
                    new_res.push_back(n_r);
                }
                m_new_chains.push_back(std::make_shared<Chain>(new_res));
            }

            BasepairOPs m_bps, m_ends;
            for(auto const & bp : m->basepairs()) {
                auto new_bp = n_ss.get_bp(n_ss.get_residue(bp->res1()->uuid()),
                                          n_ss.get_residue(bp->res2()->uuid()));
                m_bps.push_back(new_bp);
            }
            
            for(auto const & end : m->ends()) {
                int i = (int)(std::find(m->basepairs().begin(), m->basepairs().end(), end) -
                                        m->basepairs().begin());
                m_ends.push_back(m_bps[i]);
            }
            auto m_copy = std::make_shared<Motif>(m->type(), m_ends, m_new_chains);
            m_copy->basepairs(m_bps);
            m_copy->name(m->name());
            m_copy->end_ids(m->end_ids());
            if (motifs.find(m_copy->type()) == motifs.end()) {
                motifs[m_copy->type()] = MotifOPs();
            }
            motifs[m_copy->type()].push_back(m_copy);
            motifs["ALL"].push_back(m_copy);
        }
        n_ss.set_motifs(motifs);
        
        return n_ss;
    }
    
    
    inline
    MotifOPs const &
    motifs(
        String const & type) {
    
        try {
            return motifs_[type];
        }
        catch(...) {
            throw SecondaryStructureException("no motifs of type " + type + "in secondary structure");
        }
    }
    
    String
    to_str() const;
    
public: //setters
    
    inline
    void
    set_motifs(std::map<String, MotifOPs> const & motifs) { motifs_ = motifs; }
    
private:
    std::map<String, MotifOPs> motifs_;
    
};
    
typedef std::shared_ptr<SecondaryStructure> SecondaryStructureOP;
 
String
assign_end_id(
    SecondaryStructureOP const &,
    BasepairOP const &);
 

//helper functions for str_to_secondary_structure
ResidueOP
get_res_from_res_str(
    SecondaryStructure const &,
    String const &);

    
SecondaryStructure
str_to_secondary_structure(
    String const &);
    
}






#endif /* defined(__RNAMake__secondary_structure__) */
