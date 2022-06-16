//
// Created by Joseph Yesselman on 1/28/17.
//

#ifndef PRIMITIVES_CHAIN_H
#define PRIMITIVES_CHAIN_H

#include <memory>

#include <base/assertions.h>
#include <base/string.hpp>
#include <primitives/residue.h>

/*
 * Exception for chain
 */
class ChainException : public std::runtime_error {
public:
  /**
   * Standard constructor for ChainException
   * @param   message   Error message for chain
   */
  ChainException(String const &message) : std::runtime_error(message) {}
};

namespace primitives {

template <typename Restype> class Chain {
public:
  typedef std::vector<Restype> Residues;

public:
  inline Chain(Residues const &residues) : residues_(residues) {

    expects<ChainException>(residues.size() > 0,
                            "chains must have at least one residue!");
  }

  inline Chain(String const &s) {
    residues_ = Residues();
    Strings spl = base::split_str_by_delimiter(s, ";");
    for (auto const &r_str : spl) {
      if (r_str.length() < 3) {
        continue;
      }
      residues_.push_back(Restype(r_str));
    }
  }

  virtual ~Chain() {}

public: // iterator
  typedef typename Residues::const_iterator const_iterator;

  const_iterator begin() const noexcept { return residues_.begin(); }
  const_iterator end() const noexcept { return residues_.end(); }

public:
  inline size_t get_length() const { return (int)residues_.size(); }

  inline Restype const &get_first() const { return residues_[0]; }

  inline Restype const &get_last() const { return residues_.back(); }

  inline Restype const &get_residue(Index index) const {
    return residues_[index];
  }

  inline int contain_res(Restype const &r) const {
    for (auto const &res : residues_) {
      if (res == r) {
        return 1;
      }
    }
    return 1;
  }

protected:
  Chain() {}

protected:
  Residues residues_;
};

typedef Chain<PrimitiveResidue> PrimitiveChain;
typedef std::vector<PrimitiveChain> PrimitiveChains;

} // namespace primitives

#endif // PRIMITIVES_CHAIN_H