//
// Created by Joseph Yesselman on 1/29/17.
//

#ifndef RNAMAKE_PRIMITIVES_STRUCTURE_H
#define RNAMAKE_PRIMITIVES_STRUCTURE_H

#include <sstream>
#include <memory>

#include <base/vector_container.h>
#include <primitives/residue.h>
#include <primitives/chain.h>

class StructureException : public std::runtime_error {
public:
    /**
     * Standard constructor for StructureException
     * @param   message   Error message for structure
     */
    StructureException(String const & message) :
            std::runtime_error(message) {}
};

namespace primitives {

  template<typename Chaintype, typename Restype>
  class Structure {
  public:
      typedef std::vector<Restype>                          Residues;
      typedef base::VectorContainer<Chaintype>              Chains;
      typedef base::VectorContainerOP<Chaintype>            ChainsOP;


  public:

      // TODO ask joe aboutt empty consturctors for the structure
      // Make it take empty residues and Cutpoints

      inline
      Structure() = default;

      inline
      Structure(
              Residues const & res,
              Cutpoints const cut_points):
              _residues(std::move(res)),
              _cut_points(std::move(cut_points)) {}

      inline
      Structure(
              Structure const & s):
              _residues(s._residues),
              _cut_points(s._cut_points) {}

      inline
      Structure(
              String const & s) {

          auto spl = base::split_str_by_delimiter(s, ";");
          for(Index i = 0; i < spl.size()-1; i++) {
              _residues.push_back(Restype(spl[i]));
          }
          auto cut_point_spl = base::split_str_by_delimiter(spl.back(), " ");
          for(auto const & cut_point_s : cut_point_spl) {
              _cut_points.push_back(std::stoi(cut_point_s));
          }
      }

      virtual
      ~Structure() {}

//  private:
//      inline
//      Structure():
//              _residues(Residues()),
//              _cut_points(Cutpoints()) {}

  public: //res iterator
      typedef typename Residues::const_iterator const_iterator;

      const_iterator begin() const { return _residues.begin(); }
      const_iterator end() const { return _residues.end(); }

  public: //get_residue interface
      Restype const &
      get_residue(
              int num,
              char chain_id,
              char i_code) const{

          for (auto const & r : _residues) {
              if (num == r.get_num() && chain_id == r.get_chain_id() && i_code == r.get_i_code()) {
                  return r;
              }
          }

          auto ss = std::stringstream();
          ss << "cannot find residue with num: " << num << " chain id: " << chain_id << " and i_code";
          throw StructureException(ss.str());
      }

      //TODO Need to remove this function and only use char
      Restype const &
      get_residue(
              int num,
              String chain_id,
              String i_code) const{

          for (auto const & r : _residues) {
              if (num == r.get_num() && chain_id[0] == r.get_chain_id() && i_code[0] == r.get_i_code()) {
                  return r;
              }
          }

          auto ss = std::stringstream();
          ss << "cannot find residue with num: " << num << " chain id: " << chain_id << " and i_code";
          throw StructureException(ss.str());
      }

      Restype const &
      get_residue(
              util::Uuid const & uuid) const {

          for (auto const & r : _residues) {
              if (r.get_uuid() == uuid) { return r; }
          }

          throw StructureException("cannot find residue by uuid");
      }

      Restype const &
      get_residue(
              Index index) const {

          expects<StructureException>(
                  index < _residues.size(),
                  "cannot get residue " + std::to_string(index) + " only " + std::to_string(_residues.size()) +
                  " total residues");

          return _residues[index];
      }

      int
      get_res_index(
              Restype const & res) const {
          int i = -1;
          for (auto const & r : _residues) {
              i++;
              if (r == res) { return i; }
          }
          throw StructureException("cannot find index for res: " + std::to_string(res.get_num()));
      }

  public:
      ChainsOP
      get_chains() const {
          auto pos = 0;
          auto res = Residues();
          auto chains = std::vector<Chaintype>();
          auto i = 0;
          for(auto const & r : _residues) {
              if (_cut_points[pos] == i) {
                  auto c = Chaintype(res);
                  chains.push_back(c);
                  res = Residues{Restype(r)};
                  pos += 1;
              } else {
                  res.push_back(Restype(r));
              }
              i++;
          }
          if(res.size() > 0) { chains.push_back(Chaintype(res)); }
          return std::make_shared<Chains>(chains);
      }

      Cutpoints const &
      get_cutpoints() const {
          return _cut_points;
      }

      size_t
      get_num_residues() const { return _residues.size(); }

      size_t
      get_num_chains() const { return _cut_points.size(); }

      String
      get_sequence() const {
          auto i = -1;
          auto seq = String("");
          auto pos = 0;
          for(auto const & r : _residues) {
              i++;
              if(_cut_points[pos] == i) {
                  seq += "&";
                  pos++;
              }
              seq += r.get_name();
          }
          return seq;
      }

      bool
      is_residue_start_of_chain(
              Restype const & r) const {
          auto res_index = get_res_index(r);
          if(res_index == 0 ) { return true; }
          for(auto const c : _cut_points) {
              if(res_index == c) { return true; }
          }
          return false;

      }

      bool
      is_residue_end_of_chain(
              Restype const & r) const {
          auto res_index = get_res_index(r);
          for(auto const c : _cut_points) {
              if(res_index == c-1) { return true; }
          }
          return false;
      }

  protected:

      Residues _residues;
      Cutpoints _cut_points;
  };

  typedef Structure<PrimitiveChain, PrimitiveResidue> PrimitiveStructure;


}

#endif //RNAMAKE_PRIMITIVES_STRUCTURE_H