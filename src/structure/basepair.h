//
// Created by Joseph Yesselman on 12/15/17.
//

#ifndef RNAMAKE_NEW_ALL_ATOM_BASEPAIR_H
#define RNAMAKE_NEW_ALL_ATOM_BASEPAIR_H

#include <base/json.h>
#include <base/vector_container.h>
#include <math/xyz_matrix.h>
#include <math/xyz_vector.h>
#include <math/numerical.h>
#include <util/x3dna.h>
#include <primitives/basepair.h>
#include <structure/residue.h>

namespace structure {

  class Basepair : public primitives::Basepair {
  public:
      inline
      Basepair(
              util::Uuid const &res1_uuid,
              util::Uuid const &res2_uuid,
              util::Uuid const &uuid,
              primitives::BasepairType bp_type,
              base::SimpleStringCOP name,
              util::X3dnaBPType x3dna_type,
              math::Matrix const &ref_frame,
              math::Point const &center,
              math::Points const &c1_prime_coords) :
              primitives::Basepair(res1_uuid, res2_uuid, uuid, bp_type, name),
              _ref_frame(ref_frame),
              _center(center),
              _c1_prime_coords(c1_prime_coords),
              _x3dna_type(x3dna_type) {}

      inline
      Basepair(
              json::JSON &j,
              util::Uuid const &res1_uuid,
              util::Uuid const &res2_uuid,
              util::Uuid const &uuid) :
              primitives::Basepair() {
          _uuid = uuid;
          _res1_uuid = res1_uuid;
          _res2_uuid = res2_uuid;
          _center = math::Point(j[0]);
          _ref_frame = math::Matrix(j[1]);
          _c1_prime_coords = math::Points(2);
          _c1_prime_coords[0] = math::Point(j[2]);
          _c1_prime_coords[1] = math::Point(j[3]);
          _bp_type = static_cast<primitives::BasepairType>(j[4].ToInt());
          _x3dna_type = static_cast<util::X3dnaBPType>(j[5].ToInt());
          _name = std::make_shared<base::SimpleString>(j[6].ToString());
      }


  public:
      bool
      is_equal(
              Basepair const &bp,
              bool check_uuid = true) const {
          if (check_uuid) {
              if (_uuid != bp._uuid) { return false; };
              if (_res1_uuid != bp._res1_uuid) { return false; }
              if (_res2_uuid != bp._res2_uuid) { return false; }
          }

          if (!math::are_points_equal(_center, bp._center)) { return false; }
          if (!math::are_matrices_equal(_ref_frame, bp._ref_frame)) { return false; }
          if (!math::are_points_equal(_c1_prime_coords[0], bp._c1_prime_coords[0])) { return false; }
          if (!math::are_points_equal(_c1_prime_coords[1], bp._c1_prime_coords[1])) { return false; }
          return true;

      }


  public: // non const methods

      inline
      void
      move(
              math::Point const &p) {
          _center = _center + p;
          _c1_prime_coords[0] = _c1_prime_coords[0] + p;
          _c1_prime_coords[1] = _c1_prime_coords[1] + p;
      }

      inline
      void
      transform(
              math::Matrix const &r,
              math::Vector const &t,
              math::Point &dummy) {
          math::dot_vector(r, _center, dummy);
          _center = dummy + t;

          math::dot_vector(r, _c1_prime_coords[0], dummy);
          _c1_prime_coords[0] = dummy + t;

          math::dot_vector(r, _c1_prime_coords[1], dummy);
          _c1_prime_coords[1] = dummy + t;

          _ref_frame = math::dot(_ref_frame, r);
          _ref_frame.unitarize();
      }

      inline
      void
      transform(
              math::Matrix const &r,
              math::Vector const &t) {
          auto dummy = math::Point();
          transform(r, t, dummy);
      }

      inline
      void
      swap_residue_positions() {
          std::swap(_res1_uuid, _res2_uuid);
          std::swap(_c1_prime_coords[0], _c1_prime_coords[1]);
      }

      inline
      void
      invert_reference_frame() {
          _ref_frame = _ref_frame.get_flip_orientation();
      }

      inline
      void
      new_uuids(
              util::Uuid const &r1_uuid,
              util::Uuid const &r2_uuid) {
          _res1_uuid = r1_uuid;
          _res2_uuid = r2_uuid;
          _uuid = util::Uuid();
      }

  public:
      json::JSON
      get_json() const {
          return json::Array(
                  _center.get_json(), _ref_frame.get_json(), _c1_prime_coords[0].get_json(),
                  _c1_prime_coords[1].get_json(), (int) _bp_type, (int) _x3dna_type, _name->get_str());
      }


  public: // getters
      inline
      math::Matrix const &
      get_ref_frame() const { return _ref_frame; }

      inline
      math::Point const &
      get_center() const { return _center; }

      inline
      math::Points const &
      get_c1_prime_coords() const { return _c1_prime_coords; }

      inline
      math::Point const &
      get_res1_c1_prime_coord() const { return _c1_prime_coords[0]; }

      inline
      math::Point const &
      get_res2_c1_prime_coord() const { return _c1_prime_coords[1]; }

  //TODO Remove setters
  public: //setters
      inline
      void
      set_ref_frame(math::Matrix r) {
          _ref_frame = r;
      }

      inline
      void
      set_uuid(util::Uuid const & id) {
          _uuid = id;
      }

  private:
      math::Point _center;
      math::Points _c1_prime_coords;
      math::Matrix _ref_frame;
      util::X3dnaBPType _x3dna_type;

  };

  typedef std::shared_ptr<Basepair> BasepairOP;
  typedef std::vector<Basepair> Basepairs;
  typedef std::vector<BasepairOP> BasepairOPs;

  inline
  String
  generate_bp_name(
          Residue const &res1,
          Residue const &res2) {
      return primitives::generate_bp_name<Residue>(res1, res2);
  }

  primitives::BasepairType
  generate_bp_type(
          Residue const &,
          Residue const &,
          util::X3dnaBPType);

}

#endif //RNAMAKE_NEW_ALL_ATOM_BASEPAIR_H