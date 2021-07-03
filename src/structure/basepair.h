//
// Created by Joseph Yesselman on 12/15/17.
//

#ifndef RNAMAKE_NEW_ALL_ATOM_BASEPAIR_H
#define RNAMAKE_NEW_ALL_ATOM_BASEPAIR_H

#include <base/json.h>
#include <base/vector_container.h>
#include <math/xyz_matrix.h>
#include <math/xyz_vector.h>
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
              ref_frame_(ref_frame),
              center_(center),
              c1_prime_coords_(c1_prime_coords),
              x3dna_type_(x3dna_type) {}

      inline
      Basepair(
              json::JSON &j,
              util::Uuid const &res1_uuid,
              util::Uuid const &res2_uuid,
              util::Uuid const &uuid) :
              primitives::Basepair() {
          uuid_ = uuid;
          res1_uuid_ = res1_uuid;
          res2_uuid_ = res2_uuid;
          center_ = math::Point(j[0]);
          ref_frame_ = math::Matrix(j[1]);
          c1_prime_coords_ = math::Points(2);
          c1_prime_coords_[0] = math::Point(j[2]);
          c1_prime_coords_[1] = math::Point(j[3]);
          bp_type_ = static_cast<primitives::BasepairType>(j[4].ToInt());
          x3dna_type_ = static_cast<util::X3dnaBPType>(j[5].ToInt());
          name_ = std::make_shared<base::SimpleString>(j[6].ToString());
      }


  public:
      bool
      is_equal(
              Basepair const &bp,
              bool check_uuid = true) const {
          if (check_uuid) {
              if (uuid_ != bp.uuid_) { return false; };
              if (res1_uuid_ != bp.res1_uuid_) { return false; }
              if (res2_uuid_ != bp.res2_uuid_) { return false; }
          }

          if (!math::are_points_equal(center_, bp.center_)) { return false; }
          if (!math::are_matrices_equal(ref_frame_, bp.ref_frame_)) { return false; }
          if (!math::are_points_equal(c1_prime_coords_[0], bp.c1_prime_coords_[0])) { return false; }
          if (!math::are_points_equal(c1_prime_coords_[1], bp.c1_prime_coords_[1])) { return false; }
          return true;

      }


  public: // non const methods

      inline
      void
      move(
              math::Point const &p) {
          center_ = center_ + p;
          c1_prime_coords_[0] = c1_prime_coords_[0] + p;
          c1_prime_coords_[1] = c1_prime_coords_[1] + p;
      }

      inline
      void
      transform(
              math::Matrix const &r,
              math::Vector const &t,
              math::Point &dummy) {
          math::dot_vector(r, center_, dummy);
          center_ = dummy + t;

          math::dot_vector(r, c1_prime_coords_[0], dummy);
          c1_prime_coords_[0] = dummy + t;

          math::dot_vector(r, c1_prime_coords_[1], dummy);
          c1_prime_coords_[1] = dummy + t;

          ref_frame_ = math::dot(ref_frame_, r);
          ref_frame_.unitarize();
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
          std::swap(res1_uuid_, res2_uuid_);
          std::swap(c1_prime_coords_[0], c1_prime_coords_[1]);
      }

      inline
      void
      invert_reference_frame() {
          ref_frame_ = ref_frame_.get_flip_orientation();
      }

      inline
      void
      new_uuids(
              util::Uuid const &r1_uuid,
              util::Uuid const &r2_uuid) {
          res1_uuid_ = r1_uuid;
          res2_uuid_ = r2_uuid;
          uuid_ = util::Uuid();
      }

  public:
      json::JSON
      get_json() const {
          return json::Array(
                  center_.get_json(), ref_frame_.get_json(), c1_prime_coords_[0].get_json(),
                  c1_prime_coords_[1].get_json(), (int) bp_type_, (int) x3dna_type_, name_->get_str());
      }


  public: // getters
      inline
      math::Matrix const &
      get_ref_frame() const { return ref_frame_; }

      inline
      math::Point const &
      get_center() const { return center_; }

      inline
      math::Points const &
      get_c1_prime_coords() const { return c1_prime_coords_; }

      inline
      math::Point const &
      get_res1_c1_prime_coord() const { return c1_prime_coords_[0]; }

      inline
      math::Point const &
      get_res2_c1_prime_coord() const { return c1_prime_coords_[1]; }

  //TODO Remove setters
  public: //setters
      inline
      void
      set_ref_frame(math::Matrix r) {
          ref_frame_ = r;
      }

  private:
      math::Point center_;
      math::Points c1_prime_coords_;
      math::Matrix ref_frame_;
      util::X3dnaBPType x3dna_type_;

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