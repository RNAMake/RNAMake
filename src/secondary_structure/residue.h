
//
//  residue.h
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__ss_residue__
#define __RNAMake__ss_residue__

#include <stdio.h>
#include <sstream>
#include <memory>
#include <stdexcept>
#include <cassert>

#include <base/types.hpp>
#include <base/string.hpp>
#include <util/uuid.h>

#include <primitives/residue.h>
#include <secondary_structure/residue.fwd.h>

namespace secondary_structure {

  class Residue : public primitives::Residue {
  public:
      //inline
      //Residue()  {}

      inline
      Residue(
              char name,
              char dot_bracket,
              int num,
              char chain_id,
              char const & i_code,
              util::Uuid const & uuid) :
              dot_bracket_(dot_bracket),
              primitives::Residue(name, num, chain_id, i_code, uuid) {
          res_type_ = _assign_res_type(_name);
      }

      inline
      Residue(
              Residue const & r) :
              Residue(r._name, r.dot_bracket_, r._num, r._chain_id, r._i_code, r._uuid) {}

      inline
      Residue(
              String const & s) {
          Strings spl = base::split_str_by_delimiter(s, ",");
          if (spl.size() != 5) {
              throw ResidueException("invalid number of elements in string that encodes secondary_structure::Residue");
          }
          _name = spl[0][0];
          dot_bracket_ = spl[1][0];
          _num = std::stoi(spl[2]);
          _chain_id = spl[3][0];
          _i_code = spl[4][0];
          _uuid = util::Uuid();
          res_type_ = _assign_res_type(_name);
      }

      virtual
      ~Residue() {}

  public:

      inline
      bool
      operator == (
              Residue const & r) const {
          if(_name != r._name) { return false; }
          if(dot_bracket_ != r.dot_bracket_) { return false; }
          if(_num != r._num) { return false; }
          if(_chain_id != r._chain_id) { return false; }
          if(_i_code != r._i_code) { return false; }
          if(_uuid != r._uuid) { return false; }
          return true;
      }

      inline
      bool
      operator != (
              Residue const & r) const {
          return !(*this == r);
      }

  public:
      inline
      bool
      is_equal(
              Residue const & r,
              bool check_uuid = true) const {
          if(_name != r._name) { return false; }
          if(dot_bracket_ != r.dot_bracket_) { return false; }
          if(_num != r._num) { return false; }
          if(_chain_id != r._chain_id) { return false; }
          if(_i_code != r._i_code) { return false; }
          if(check_uuid && _uuid != r._uuid) { return false; }
          return true;
      }


  public:
      inline
      String
      get_str() const {
          std::stringstream ss;
          ss << _name << "," << dot_bracket_ << "," << _num << "," << _chain_id << "," << _i_code;
          return ss.str();
      }

  public: //getters
      inline
      char
      get_dot_bracket() const { return dot_bracket_; }

      inline
      int
      get_res_type() const { return res_type_; }

  public: //setters

      inline
      void
      set_name(char name) {
          _name = name;
          res_type_ = _assign_res_type(_name);
      }

  private:
      int
      _assign_res_type(char name) {
          if (name == 'A') { return 0; }
          else if (name == 'C') { return 1; }
          else if (name == 'G') { return 2; }
          else if (name == 'U') { return 3; }
          else if (name == 'T') { return 3; }
          else if (name == 'N') { return -1; }
          else {
              throw ResidueException(
                      "in sstruct::Residue encountered a unknown name: " + std::string(1, name));
          }
      }

  private:
      //A=0,C=1,G=2,U=3,T=3
      int res_type_;
      char dot_bracket_;

  };

/*struct res_less_than_key {
    inline
    bool
    operator() (ResidueOP const & r1, ResidueOP const & r2) {
        return (r1->num() < r2->num());
    }
};*/

  typedef std::vector<Residue> Residues;

} //secondary_structure

#endif /* defined(__RNAMake__ss_residue__) */