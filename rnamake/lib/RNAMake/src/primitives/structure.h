//
// Created by Joseph Yesselman on 1/29/17.
//

#ifndef RNAMAKE_PRIMITIVES_STRUCTURE_H
#define RNAMAKE_PRIMITIVES_STRUCTURE_H

#include "primitives/residue.h"
#include "primitives/chain.h"

class StructureException : public std::runtime_error {
public:
    /**
     * Standard constructor for StructureException
     * @param   message   Error message for structure
     */
    StructureException(String const & message) :
            std::runtime_error(message) {}
};

namespace primitives {

template<typename Chaintype, typename Restype>
class Structure {
public:
    typedef std::shared_ptr<Restype> ResidueOP;
    typedef std::vector<ResidueOP>   ResidueOPs;

    typedef std::shared_ptr<Chaintype> ChainOP;
    typedef std::vector<ChainOP>       ChainOPs;

public:
    inline
    Structure(ChainOPs const & chains):
            chains_(chains),
            residues_(ResidueOPs()) {
        for(auto const & c : chains_) {
            for(auto const & r : *c) { residues_.push_back(r); }
        }
    }

    virtual
    ~Structure() {}

public: //res iterator
    typedef typename ResidueOPs::iterator iterator;
    typedef typename ResidueOPs::const_iterator const_iterator;

    iterator begin() { return residues_.begin(); }
    iterator end()   { return residues_.end(); }

    const_iterator begin() const { return residues_.begin(); }
    const_iterator end()   const { return residues_.end(); }

public: //chain iterator
    typedef typename ChainOPs::iterator chain_iterator;
    typedef typename ChainOPs::const_iterator chain_const_iterator;

    chain_iterator chain_begin() { return chains_.begin(); }
    chain_iterator chain_end()   { return chains_.end(); }

    chain_const_iterator chain_begin() const { return chains_.begin(); }
    chain_const_iterator chain_end()   const { return chains_.end(); }

public:
    ResidueOP const
    get_residue(
            int const & num,
            String const & chain_id,
            String const & i_code) {

        for (auto const & r : residues_) {
            if (num == r->num() && chain_id == *(r->chain_id()) && i_code == *(r->i_code())) {
                return r;
            }
        }
        return ResidueOP(nullptr);
    }

    ResidueOP const
    get_residue(
            Uuid const & uuid) {

        for (auto const & r : residues_) {
            if (r->uuid() == uuid) { return r; }
        }

        return ResidueOP(nullptr);
    }

    ResidueOP const &
    get_residue(
            int index) {
        return residues_[index];
    }

    int
    get_res_index(ResidueOP const & res) {
        int i = -1;
        for(auto const & r : residues_) {
            i++;
            if(r == res) { return i; }
        }
        throw StructureException("cannot find index for res: " + std::to_string(res->num()));
    }


public:
    ChainOP const &
    get_chain(int i) { return chains_[i]; }

    size_t
    num_residues() { return residues_.size(); }

    size_t
    num_chains() { return chains_.size(); }


protected:
    Structure() {}

protected:

    ChainOPs chains_;
    ResidueOPs residues_;
};

}

#endif //RNAMAKE_PRIMITIVES_STRUCTURE_H
