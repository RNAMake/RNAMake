//
// Created by Joseph Yesselman on 1/26/17.
//

#ifndef PRIMITIVES_RESIDUE_H
#define PRIMITIVES_RESIDUE_H

#include <base/types.h>
#include <util/uuid.h>

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
        /**
         * constructor for Residue class
         * @param   name        residue name (A, G, C, T)
         * @param   num         residue num
         * @param   chain_id    what chain does this residue belong to ("A", "B")
         * @param   uuid        residue unique indentifier
         */
        inline
        Residue(
                char name,
                int num,
                char chain_id,
                char i_code,
                util::Uuid const & uuid):
                name_(name),
                num_(num),
                chain_id_(chain_id),
                i_code_(i_code),
                uuid_(uuid) {}

        /**
          * copy construtor for residue class
          * @param  r   residue to be copied from
          */
        inline
        Residue(
                Residue const & r):
                name_(r.name_),
                num_(r.num_),
                chain_id_(r.chain_id_),
                uuid_(r.uuid_),
                i_code_(r.i_code_) {}

        virtual
        ~Residue() {}

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
        operator==(Residue const & other) const {
            return uuid_ == other.uuid_;
        }

        inline
        bool
        operator!=(Residue const & other) const {
            return uuid_ != other.uuid_;
        }

    public: //getters

        /**
         * getter the chain_id, i.e. "A", "B", the id of the chain this residue belongs to
         */
        inline
        char
        get_chain_id() const { return chain_id_; }

        /**
         * getter for the name of the residue, i.e. "A", "G" etc
         */
        inline
        char
        get_name() const { return name_; }

        /**
         * getter for the residue num
         */
        inline
        int
        get_num() const { return num_; }

        /**
         * getter for the residue insertion code
         */
        inline
        char
        get_i_code() const { return i_code_; }

        /**
        * getter for residue unique indentifier
        */
        inline
        util::Uuid const &
        get_uuid() const { return uuid_; }

        virtual
        String
        get_str() const { return String(""); }

    protected:

        char name_;
        int num_;
        char chain_id_;
        char i_code_;
        util::Uuid uuid_;

    };

    typedef Residue                          PrimitiveResidue;
    typedef std::vector<Residue>             PrimitiveResidues;

}
#endif //PRIMITIVES_RESIDUE_H