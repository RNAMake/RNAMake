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
  inline Structure() : _residues(Residues()), _cut_points(Cutpoints()) {}

  inline Structure(Residues &res, Cutpoints &cut_points)
      : _residues(std::move(res)), _cut_points(std::move(cut_points)) {}

  ~Structure() = default;

public: // res iterator
  typedef typename Residues::const_iterator const_iterator;

  const_iterator begin() const { return _residues.begin(); }
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

  Residue const &get_residue(util::Uuid const &uuid) const {
    for (auto const &r : _residues) {
      if (r.get_uuid() == uuid) {
        return r;
      }
    }
    throw StructureException("cannot find residue by uuid");
  }

  Residue const &get_residue(Index index) const {
    /*expects<StructureException>(
        index < _residues.size(),
        "cannot get residue " + std::to_string(index) + " only " +
            std::to_string(_residues.size()) + " total residues");   */

    return _residues[index];
  }

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

  [[nodiscard]] const Cutpoints &get_cutpoints() const { return _cut_points; }

  [[nodiscard]] size_t get_num_residues() const { return _residues.size(); }

  [[nodiscard]] size_t get_num_chains() const { return _cut_points.size(); }

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
  void move(const math::Vector3 &p) {
    for (auto &r : _residues) {
      r.move(p);
    }
  }

  void transform(const math::RotandTrans &rt) {
    for (auto &r : _residues) {
      r.transform(rt);
    }
  }

private:
  Residues _residues;
  Cutpoints _cut_points;
};

}

#endif // RNAMAKE_SRC_STRUCTURE_BASE_STRUCTURE_HPP_
