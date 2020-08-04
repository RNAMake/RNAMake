//
// Created by Joseph Yesselman on 2019-04-24.
//

#ifndef RNAMAKE_NEW_STRUCTURE_BEADS_H
#define RNAMAKE_NEW_STRUCTURE_BEADS_H

#include "math/xyz_vector.h"

namespace structure {

/**
 * BeadType is an ENUM type. This is to specify which center of atoms each bead
 * represents.
 *
 * Phosphate (0):  P, OP1, OP2\n
 * Sugar (1):  O5',C5',C4',O4',C3',O3',C1',C2',O2'\n
 * Base  (2):  All remaining atoms
 */
enum class BeadType {
    PHOS = 0,
    SUGAR = 1,
    BASE = 2
};

/**
 * Bead class stores information related to keeping track of steric clashes
 * between residues during building. They are never used outside the Residue class
 *
 */
class Bead {
public:
    // TODO I do not like these empty constructors -- JDY
    /**
     * Empty constructor for Bead object.
     */
    Bead() :
            _center(math::Point(0, 0, 0)),
            _btype(BeadType::PHOS) {}

    /**
     * Standard constructor for Bead object.
     * @param   btype   type of bead (PHOS, SUGAR or BASE)
     * @param   center  the average 3D position of all atoms represented by bead
     */
    inline
    Bead(
            math::Point const & center,
            BeadType const btype) :
            _center(center),
            _btype(btype) {}


    inline
    Bead(
            String const & s) {
        auto spl = base::split_str_by_delimiter(s, ",");
        _center = math::vector_from_str(spl[0]);
        _btype = BeadType(std::stoi(spl[1]));
    }

    /**
     * Copy constructor
     * @param   b   Bead object copying from
     */
    inline
    Bead(
            Bead const & b) :
            _center(b._center),
            _btype(b._btype) {}

public:

    inline
    double
    distance(
            Bead const & b) const {
        return b._center.distance(_center);
    }

public: //accessors

    /**
     * Accessor for _center
     */
    inline
    math::Point
    center() const {
        return _center;
    }

    /**
     * Accessor for _btype
     */
    inline
    BeadType
    btype() const {
        return _btype;
    }

private:
    /**
     * private variable for the 3D coordinates of the center of atoms the bead represents
     */
    math::Point _center;

    /**
     * private variable of the type of the bead PHOS, SUGAR or BASE)
     */
    BeadType _btype;

};

typedef std::vector<Bead> Beads;

}

#endif //RNAMAKE_NEW_STRUCTURE_BEADS_H
