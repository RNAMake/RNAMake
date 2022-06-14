//
// Created by Joseph Yesselman on 11/12/17.
//

#ifndef RNAMAKE_NEW_BEAD_H
#define RNAMAKE_NEW_BEAD_H

#include <math/vector_3.hpp>
#include <math/matrix_3x3.hpp>

namespace util {

/**
 * BeadType is an ENUM type. This is to specify which center of atoms each bead
 * represents.
 *
 * Phosphate (0):  P, OP1, OP2\n
 * Sugar (1):  O5',C5',C4',O4',C3',O3',C1',C2',O2'\n
 * Base  (2):  All remaining atoms
 */
    enum class BeadType {
        PHOS,
        SUGAR,
        BASE,
        CALPHA,
        MCENTER // center of small molecule
    };

/**
 * Exception for beads
 */
    class BeadException : public std::runtime_error {
    public:
        /**
         * Standard constructor for BeadException
         * @param   message   Error message for bead
         */
        BeadException(String const & message) :
                std::runtime_error(message) {}
    };


/**
 * Bead class stores information related to keeping track of steric clashes
 * between residues during building. They are never used outside the Residue class
 *
 */
    class Bead {
    public:
        /**
         * Standard constructor for Bead object.
         * @param   btype   type of bead (PHOS, SUGAR or BASE)
         * @param   center  the average 3D position of all atoms represented by bead
         */
        inline
        Bead(
                math::Vector3 const & center,
                BeadType const bead_type):
                center_ ( center ),
                bead_type_ ( bead_type )  {}

        inline
        Bead(
                String const & s) {
            auto spl     = base::split_str_by_delimiter(s, ",");
            center_      = math::Vector3(spl[0]);
            bead_type_   = BeadType(std::stoi(spl[1]));
        }

        /**
         * Copy constructor
         * @param   b   Bead object copying from
         */
        inline
        Bead(
                Bead const & b):
                center_(b.center_),
                bead_type_(b.bead_type_) {}

    public:

        inline
        double
        distance(
                Bead const & b) const {
            return b.center_.distance(center_);
        }

        inline
        void
        move(
                math::Vector3 const & p) {
            center_ = center_ + p;
        }

        inline
        void
        transform(
                math::Matrix3x3 const & r,
                math::Vector3 const & t,
                math::Vector3 & dummy) {
            math::dot_vector(r, center_, dummy);
            center_ = dummy + t;
        }

        inline
        void
        transform(
                math::Matrix3x3 const & r,
                math::Vector3 const & t) {
            auto dummy = math::dot_vector(r, center_);
            center_ = dummy + t;
        }

    public: //getters

        /**
         * Accessor for center_
         */
        inline
        math::Vector3 const &
        get_center() const { return center_; }

        /**
         * Accessor for btype_
         */
        inline
        BeadType
        get_type() const { return bead_type_; }

        String
        get_str() const {
            return center_.get_str() + "," + std::to_string((int)bead_type_);
        }

        String
        get_type_name() {
            if     (bead_type_ == BeadType::PHOS)  { return "PHOSPHATE"; }
            else if(bead_type_ == BeadType::SUGAR) { return "SUGAR"; }
            else if(bead_type_ == BeadType::BASE)  { return "BASE"; }
            else { throw BeadException("unknown bead type"); }
        }


    private:
        /**
         * private variable for the 3D coordinates of the center of atoms the bead represents
         */
        math::Vector3 center_;

        /**
         * private variable of the type of the bead PHOS, SUGAR or BASE)
         */
        BeadType bead_type_;

    };

    typedef std::vector<Bead> Beads;

}


#endif //RNAMAKE_NEW_BEAD_H