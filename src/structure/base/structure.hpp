//
// Created by Joe Yesselman on 6/30/22.
//

#ifndef RNAMAKE_SRC_STRUCTURE_BASE_STRUCTURE_HPP_
#define RNAMAKE_SRC_STRUCTURE_BASE_STRUCTURE_HPP_

#include <structure/base/chain.hpp>

namespace structure::base {

template <typename Chain, typename Residue> class Structure {
public:
  typedef std::vector<Residue> Residues;
  typedef std::vector<Chain> Chains;
  typedef std::vector<Chain> ChainsOP;

public:
  /// @brief - constructor
  inline Structure() : _residues(Residues()), _cut_points(Cutpoints()) {}

  inline Structure(Residues &res, Cutpoints &cut_points)
      : _residues(std::move(res)), _cut_points(std::move(cut_points)) {}

  /// @brief - deconstructor
  ~Structure() = default;

public: // res iterator
  typedef typename Residues::const_iterator const_iterator;

  /// @brief - gets the start of the residue chain
  const_iterator begin() const { return _residues.begin(); }

  /// @brief - gets the end of the residue chain
  const_iterator end() const { return _residues.end(); }

public: // get_residue interface
  // TODO Need to remove this function and only use char
  Residue const &get_residue(int num, const String &chain_id,
                             char i_code) const {
    for (auto const &r : _residues) {
      if (num == r.get_num() && chain_id == r.get_chain_id() &&
          i_code == r.get_i_code()) {
        return r;
      }
    }
    throw StructureException("cannot find residue!");

    /*auto ss = std::stringstream();
    ss << "cannot find residue with num: " << num << " chain id: " << chain_id
       << " and i_code";
    throw StructureException(ss.str());   */
  }

  /// @brief - gets a residue by the UUID
  Residue const &get_residue(util::Uuid const &uuid) const {
    for (auto const &r : _residues) {
      if (r.get_uuid() == uuid) {
        return r;
      }
    }
    throw StructureException("cannot find residue by uuid");
  }

  /// @brief - gets a residue by the index
  Residue const &get_residue(Index index) const {
    /*expects<StructureException>(
        index < _residues.size(),
        "cannot get residue " + std::to_string(index) + " only " +
            std::to_string(_residues.size()) + " total residues");   */

    return _residues[index];
  }

  /// @brief - gets the index of a residue
  int get_res_index(Residue const &res) const {
    int i = -1;
    for (auto const &r : _residues) {
      i++;
      if (r == res) {
        return i;
      }
    }
    throw StructureException("cannot find index for res: " +
                             std::to_string(res.get_num()));
  }

public:
  /// @brief -
  ChainsOP get_chains() const {
    auto pos = 0;
    auto res = Residues();
    auto chains = std::vector<Chain>();
    auto i = 0;
    for (auto const &r : _residues) {
      if (_cut_points[pos] == i) {
        auto c = Chain(res);
        chains.push_back(c);
        res = Residues{Residue(r)};
        pos += 1;
      } else {
        res.push_back(Residue(r));
      }
      i++;
    }
    if (res.size() > 0) {
      chains.push_back(Chain(res));
    }
    return std::make_shared<Chains>(chains);
  }

  /// @brief - retuns cutpoints
  [[nodiscard]] const Cutpoints &get_cutpoints() const { return _cut_points; }

  /// @brief - counts and returns the number of residues
  [[nodiscard]] size_t get_num_residues() const { return _residues.size(); }

  /// @brief - counts and returns the number of chains
  [[nodiscard]] size_t get_num_chains() const { return _cut_points.size(); }

  /// @brief -
  [[nodiscard]] String get_sequence() const {
    auto i = -1;
    auto seq = String("");
    auto pos = 0;
    for (auto const &r : _residues) {
      i++;
      if (_cut_points[pos] == i) {
        seq += "&";
        pos++;
      }
      seq += r.get_name();
    }
    return seq;
  }

  /// @brief - checks if the residue is the start of a chain
  bool is_residue_start_of_chain(Residue const &r) const {
    auto res_index = get_res_index(r);
    if (res_index == 0) {
      return true;
    }
    for (auto const c : _cut_points) {
      if (res_index == c) {
        return true;
      }
    }
    return false;
  }

  /// @brief - checks if the residue is the end of a chain
  bool is_residue_end_of_chain(Residue const &r) const {
    auto res_index = get_res_index(r);
    for (auto const c : _cut_points) {
      if (res_index == c - 1) {
        return true;
      }
    }
    return false;
  }

public:
  /// @brief - moves each residue in a chain by the specified position vector
  void move(const math::Vector3 &p) {
    for (auto &r : _residues) {
      r.move(p);
    }
  }

  /// @brief - rotates each residue in a chain by the specified rotation matrix
  void rotate(const math::Matrix3x3 &rot) {
    for (auto &r : _residues) {
      r.rotate(rot);
    }
  }

private:
  Residues _residues;
  Cutpoints _cut_points;
};

} // namespace structure::base

#endif // RNAMAKE_SRC_STRUCTURE_BASE_STRUCTURE_HPP_
