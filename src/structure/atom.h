#ifndef __RNAMake__atom__
#define __RNAMake__atom__

#include <stdio.h>

//RNAMake Headers
#include <base/types.h>
#include <base/simple_string.h>
#include <base/json.h>
#include <math/numerical.h>
#include <math/xyz_vector.h>
#include <math/xyz_matrix.h>


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
namespace structure {
  class Atom {
  public:

      /**
       * Standard constructor for Atom object.
       * @param   name    name of atom
       * @param   coords  3d coordinates of atom's position
       */
      inline
      Atom(
              String const &name,
              math::Point const &coords) :
              name_(name),
              coords_(coords) {}

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
              String const &s) {

          auto spl = base::split_str_by_delimiter(s, " ");
          name_ = spl[0];
          coords_ = math::Point(std::stof(spl[1]), std::stof(spl[2]), std::stof(spl[3]));
      }


//        inline
//        Atom(
//                json::JSON & j) {
//            name_ = std::make_shared<base::SimpleString const>(j[0].ToString());
//            coords_ = math::Point(j[1]);
//        }

      /**
       * Copy constructor
       * @param   a   atom object to from
       */
      inline
      Atom(
              Atom const &a) :
              name_(a.name_),
              coords_(a.coords_) {}

  public:
      inline
      bool
      operator==(
              Atom const &a) const {
          if (name_ != a.name_) { return false; }
          if (!math::are_points_equal(coords_, a.coords_)) { return false; }
          return true;
      }

      inline
      bool
      operator!=(
              Atom const &a) const {
          return !(*this == a);
      }


  public:

      /**
       * Strigifies atom object
       * @code
       *  auto a = Atom("P", Point(0, 1, 2));
       *  std::cout << a.to_str() << std::endl;
       *  //EXPECTED OUTPUT
       *  "H1 0.0 1.0 2.0"
       * @endcode
       */
      String
      get_str() const;

      json::JSON
      get_json() const;

      /**
       * Strigifies atom into PDB format
       * @param   acount  the number of the atom, default=1
       *
       * @code
       *  auto a = Atom("P", Point(0, 1, 2));
       *  std::cout << a.to_pdb_str() << std::endl;
       *  //EXPECTED OUTPUT
       *  "ATOM      1  P   C   A   1       1.000   2.000   3.000  1.00 62.18           P
       * @endcode
       */
      String
      to_pdb_str(int) const;

      /**
       * @param p xyz coords to move atom by
       *
       * @code
       * @endcode
       */
      inline
      void
      move(math::Point const &p) {
          coords_ = coords_ + p;
      }

      inline
      void
      transform(
              math::Matrix const &r,
              math::Vector const &t,
              math::Point &dummy) {
          math::dot_vector(r, coords_, dummy);
          coords_ = dummy + t;
      }

      inline
      void
      rename(String const &name) { name_ = name; }

      //TODO Temoporary needs to be removed!
      inline
      void
      set_coords(math::Point const &coords) { coords_ = coords; }

      //        inline
      //        void
      //        transform(
      //                math::Matrix const & r,
      //                math::Vector const & t) {
      //            auto dummy = math::dot_vector(r, coords_);
      //            coords_ = dummy + t;
      //        }

  public: //accessors

      /**
       * Accessor for name_
       */
      inline
      String
      get_name() const { return name_; }

//         inline
//         String
//         get_str_name() const { return name_->get_str(); }

      /**
       * Accessor for coords_
       */
      inline
      math::Point const &
      get_coords() const { return coords_; }


        inline
        double
        get_x() const { return coords_.get_x(); }

        inline
        double
        get_y() const { return coords_.get_y(); }

        inline
        double
        get_z() const { return coords_.get_z(); }

  private:
      /**
       * private variable of name of atom
       */
      //TODO Figure out a way to stop string copying
      String name_;

      /**
       * private variable of 3D coordinates of atom
       */
      math::Point coords_;

  };

/**
 * Shared pointer typedef for Atom. Only use shared pointers!
 */
  typedef std::shared_ptr<Atom> AtomOP;

/**
 * Typedef of a vector of shared pointer atoms, only used this.
 */
  typedef std::vector<AtomOP> AtomOPs;

/**
 * Typedef for multiple Atoms
 */
  typedef std::vector<Atom> Atoms;

}


#endif /* defined(__RNAMake__atom__) */