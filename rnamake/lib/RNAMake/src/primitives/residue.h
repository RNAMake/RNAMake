//
// Created by Joseph Yesselman on 1/26/17.
//

#ifndef PRIMITIVES_RESIDUE_H
#define PRIMITIVES_RESIDUE_H

#include "base/types.h"
#include "util/uuid.h"

/*
 * Exception for residues
 */
class ResidueException : public std::runtime_error {
public:
    /**
     * Standard constructor for ResidueException
     * @param   message   Error message for residue
     */
    ResidueException(String const & message) :
            std::runtime_error(message) {}
};

namespace primitives {

class Residue {
public:
    inline
    Residue(
            char name,
            int num,
            char chain_id):
            name_(name),
            num_(num),
            chain_id_(chain_id),
            i_code_(' '),
            uuid_(Uuid()) {}

    inline
    Residue(
            char name,
            int num,
            char chain_id,
            char i_code):
            name_(name),
            num_(num),
            chain_id_(chain_id),
            i_code_(i_code),
            uuid_(Uuid()) {}
    inline
    Residue(
            char name,
            int num,
            char chain_id,
            char i_code,
            Uuid const & uuid):
            name_(name),
            num_(num),
            chain_id_(chain_id),
            i_code_(i_code),
            uuid_(uuid) {}

protected:
    // let derived class setup members
    Residue() {}

public:

    /**
    * equal operator checks whether the unique indentifier is the same
    * @param   r   another residue to check if its the same
    */
    inline
    bool
    operator==(Residue & other) {
        return uuid_ == other.uuid_;
    }

    inline
    bool
    operator!=(Residue & other) const {
        return uuid_ != other.uuid_;
    }

public: //getters

    /**
    * getter for the name of the residue, i.e. "A", "G" etc
    */
    inline
    char
    name() { return name_; }

    /**
    * getter for the residue num
    */
    inline
    int
    num() { return num_; }

    /**
     * getter the chain_id, i.e. "A", "B", the id of the chain this residue belongs to
     */
    inline
    char
    chain_id() { return chain_id_; }

    /**
    * getter for the residue insertion code
    */
    inline
    char
    i_code() { return i_code_; }

    /**
    * getter for residue unique indentifier
    */
    inline
    Uuid const &
    uuid() const { return uuid_; }

protected:

    char name_;
    int num_;
    char chain_id_;
    char i_code_;
    Uuid uuid_;

};

typedef std::shared_ptr<Residue> ResidueOP;
typedef std::vector<ResidueOP>   ResidueOPs;

}
#endif //PRIMITIVES_RESIDUE_H
