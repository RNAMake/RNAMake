
#ifndef __RNAMake__atom__
#define __RNAMake__atom__

#include <stdio.h>

//RNAMake Headers
#include "base/types.h"
#include "math/xyz_vector.h"
#include "math/transform.h"

namespace structure {

/**
 * Stores atomic information from pdb file, design is to be extremely
 * lightweight only storing the atom name and coordinates.
 *
 * Example Usage:
 *
 * @code
 *  // creation
 *  auto a = Atom("P", Point(0, 1, 2));
 *
 *  //copy
 *  auto a2 = Atom(a);
 * @endcode
 */

class Atom {
public: // constructors

    /**
     * Standard constructor for Atom object.
     * @param   name    name of atom
     * @param   coords  3d coordinates of atom's position
     */
    inline
    Atom(
            String const & name,
            math::Point const & coords):
            _name(name),
            _coords(coords) {}

    /**
     * Construction from String, used in reading data from files
     * @param   s   string generated from to_str()
     * @see to_str()
     * 
     * Example Usage:
     * @code
     *  auto a = Atom("P", Point(0, 1, 2));
     *  auto s = a.to_str();
     *  auto a2 = Atom(s);
     * @endcode
     */
    inline
    Atom(
            String const & s) {

        auto spl = base::split_str_by_delimiter(s, " ");
        _name = spl[0];
        _coords = math::Point(std::stof(spl[1]), std::stof(spl[2]), std::stof(spl[3]));
    }

    /**
     * Copy constructor
     * @param   a   atom object to from
     */
    inline
    Atom(
            Atom const & a):
            _name(a._name),
            _coords(a._coords) {}

public: // output functions

    /**
     * Strigifies atom object
     * @code
     *  auto a = Atom("P", math::Point(0, 1, 2));
     *  std::cout << a.to_str() << std::endl;
     *  //EXPECTED OUTPUT
     *  "H1 0.0 1.0 2.0"
     * @endcode
     */
    String
    to_str();

    /**
     * Strigifies atom into PDB format
     * @param   acount  the number of the atom, default=1
     *
     * @code
     *  auto a = Atom("P", math::Point(0, 1, 2));
     *  std::cout << a.to_pdb_str() << std::endl;
     *  //EXPECTED OUTPUT
     *  "ATOM      1  P   C   A   1       1.000   2.000   3.000  1.00 62.18           P
     * @endcode
     */
    String
    to_pdb_str(int);

public: // non const methods

    inline
    void
    move(
            math::Point const & p) {
        _coords = _coords + p;
    }

    inline
    void
    transform(
            math::Transform const & t) {
        math::Matrix r = t.rotation().transpose();
        math::Point trans = t.translation();
        auto new_coords = math::Point();
        math::dot_vector(r, _coords, new_coords);
        new_coords += trans;
        _coords = new_coords;
    }

public: //accessors

    /**
     * Accessor for _name
     */
    inline
    String const &
    name() const {
        return _name;
    }

    /**
     * Accessor for _coords
     */
    inline
    math::Point const
    coords() const {
        return _coords;
    }


public: // setters

    /**
     * Setter for _coords
     */
    inline
    void
    coords(
            math::Point const & ncoords) {
        _coords = ncoords;
    }

    /**
     * Setter for _name
     */
    inline
    void
    name(
            String const & nname) {
        _name = nname;
    }

private:
    /**
     * private variable of name of atom
     */
    String _name;

    /**
     * private variable of 3D coordinates of atom
     */
    math::Point _coords;

};

/**
 * Shared pointer typedef for Atom. Only use shared pointers!
 */
typedef std::shared_ptr<Atom> AtomOP;

/**
 * Typedef of a vector of shared pointer atoms, only used this.
 */
typedef std::vector<AtomOP> AtomOPs;

}


#endif /* defined(__RNAMake__atom__) */
