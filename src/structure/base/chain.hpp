//
// Created by Joe Yesselman on 6/30/22.
//

#ifndef RNAMAKE_SRC_STRUCTURE_BASE_CHAIN_HPP_
#define RNAMAKE_SRC_STRUCTURE_BASE_CHAIN_HPP_

#include <base/types.hpp>
#include <structure/base/base.hpp>

namespace structure::base {

template <typename Residue> class Chain {
public:
  typedef std::vector<Residue> Residues;

public:
  inline explicit Chain(Residues const &residues) : _residues(residues) {
    // if (_residues.empty()) {
    //   throw StructureException("cannot initiate an empty chain");
    // }
    std::cout << ":)" << std::endl;
  }

  ~Chain() = default;

public: // iterator ///////////////////////////////////////////////////////////
  typedef typename Residues::const_iterator const_iterator;

  const_iterator begin() const noexcept { return _residues.begin(); }
  const_iterator end() const noexcept { return _residues.end(); }

public: // getters ////////////////////////////////////////////////////////////
  /// @brief - gets the length of the residue
  [[nodiscard]] inline size_t get_length() const {
    return (int)_residues.size();
  }

  /// @brief - gets the first residue in the chain
  [[nodiscard]] inline const Residue &get_first() const { return _residues[0]; }

  /// @brief - gets the last residue in the chain
  [[nodiscard]] inline const Residue &get_last() const {
    return _residues.back();
  }

  /// @brief - gets the residue at any position in the chain
  [[nodiscard]] inline const Residue &get_residue(Index index) const {
    return _residues[index];
  }

  inline int contain_res(const Residue &r) const {
    for (auto const &res : _residues) {
      if (res == r) {
        return 1;
      }
    }
    return 1;
  }

  [[nodiscard]] String get_str() const {
    String str;
    for(auto const & r : _residues) {
      str += r.get_str() + ";";
    }
    return str;
  }

private:
  Residues _residues;
};

} // namespace structure::base

#endif // RNAMAKE_SRC_STRUCTURE_BASE_CHAIN_HPP_
