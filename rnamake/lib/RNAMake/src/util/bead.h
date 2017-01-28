//
// Created by Joseph Yesselman on 1/26/17.
//

#ifndef RNAMAKE_BEAD_H
#define RNAMAKE_BEAD_H

#include "math/xyz_vector.h"
#include "math/xyz_matrix.h"
#include "math/transform.h"

/**
 * BeadType is an ENUM type. This is to specify which center of atoms each bead
 * represents.
 *
 * Phosphate (0):  P, OP1, OP2\n
 * Sugar (1):  O5',C5',C4',O4',C3',O3',C1',C2',O2'\n
 * Base  (2):  All remaining atoms
 */
enum BeadType {
    PHOS = 0, SUGAR = 1, BASE = 2
};

/**
 * Bead class stores information related to keeping track of steric clashes
 * between residues during building. They are never used outside the Residue class
 *
 */
class Bead {
public:
    /**
     * Empty constructor for Bead object.
     */
    Bead() :
            center_(Point(0, 0, 0)),
            btype_(BeadType(0)) {}

    /**
     * Standard constructor for Bead object.
     * @param   btype   type of bead (PHOS, SUGAR or BASE)
     * @param   center  the average 3D position of all atoms represented by bead
     */
    inline
    Bead(
            Point const & center,
            BeadType const btype):
            center_(center),
            btype_(btype) {}


    inline
    Bead(
            String const & s) {
        auto spl = split_str_by_delimiter(s, ",");
        center_ = vector_from_str(spl[0]);
        btype_ = BeadType(std::stoi(spl[1]));
    }

    /**
     * Copy constructor
     * @param   b   Bead object copying from
     */
    inline
    Bead(
            Bead const & b) :
            center_(b.center_),
            btype_(b.btype_) {}

public:

    inline
    double
    distance(Bead const & b) const {
        return b.center_.distance(center_);
    }

    inline
    void
    move(Point const & p) {
        center_ = center_ + p;
    }

    inline
    void
    transform(Transform const & t) {
        dot_vector(t.rotation().transpose(), center_);
        center_ = center_ + t.translation();
    }

    inline
    void
    fast_transform(
            Matrix const & r,
            Vector const & t) {
        dot_vector(r, center_);
        center_ = center_ + t;
    }

    String
    to_str() const;

public: //accessors

    /**
     * Accessor for center_
     */
    inline
    Point
    center() const { return center_; }

    /**
     * Accessor for btype_
     */
    inline
    BeadType
    btype() const { return btype_; }

private:
    /**
     * private variable for the 3D coordinates of the center of atoms the bead represents
     */
    Point center_;

    /**
     * private variable of the type of the bead PHOS, SUGAR or BASE)
     */
    BeadType btype_;

};

typedef std::vector<Bead> Beads;


#endif //RNAMAKE_BEAD_H
