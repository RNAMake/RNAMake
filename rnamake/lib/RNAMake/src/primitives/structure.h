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
    typedef std::vector<ResidueOP> ResidueOPs;

    typedef std::shared_ptr<Chaintype> ChainOP;
    typedef std::vector<ChainOP> ChainOPs;

public:
    inline
    Structure(
            ResidueOPs const & res,
            Ints const chain_cuts):
            residues_(res),
            chain_cuts_(chain_cuts) {}

    virtual
    ~Structure() {}

public: //res iterator
    typedef typename ResidueOPs::iterator iterator;
    typedef typename ResidueOPs::const_iterator const_iterator;

    iterator begin() { return residues_.begin(); }
    iterator end() { return residues_.end(); }

    const_iterator begin() const { return residues_.begin(); }
    const_iterator end() const { return residues_.end(); }

public: //get_residue interface
    ResidueOP const
    get_residue(
            int num,
            char chain_id,
            char i_code) const{

        for (auto const & r : residues_) {
            if (num == r->num() && chain_id == r->chain_id() && i_code == r->i_code()) {
                return r;
            }
        }
        return ResidueOP(nullptr);
    }

    ResidueOP const
    get_residue(
            Uuid const & uuid) const {

        for (auto const & r : residues_) {
            if (r->uuid() == uuid) { return r; }
        }

        return ResidueOP(nullptr);
    }

    ResidueOP const &
    get_residue(
            int index) const {
        return residues_[index];
    }

    int
    get_res_index(ResidueOP const & res) {
        int i = -1;
        for (auto const & r : residues_) {
            i++;
            if (r == res) { return i; }
        }
        throw StructureException("cannot find index for res: " + std::to_string(res->num()));
    }

public:
    size_t
    num_residues() { return residues_.size(); }

    size_t
    num_chains() { return chain_cuts_.size(); }

    String
    sequence() {
        auto i = -1;
        auto seq = String("");
        auto pos = 0;
        for(auto const & r : residues_) {
            i++;
            if(chain_cuts_[pos] == i) {
                seq += "&";
                pos++;
            }
            seq += r->name();
        }
        return seq;
    }


protected:
    Structure() {}

protected:

    ResidueOPs residues_;
    Ints chain_cuts_;
};

}

#endif //RNAMAKE_PRIMITIVES_STRUCTURE_H
