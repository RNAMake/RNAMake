//
// Created by Joseph Yesselman on 3/22/19.
//

#ifndef RNAMAKE_NEW_SOLUTION_FILTER_H
#define RNAMAKE_NEW_SOLUTION_FILTER_H

#include <memory>

#include <base/string.hpp>
#include <base/types.hpp>

namespace motif_search {

class SolutionFilter {
public:
  SolutionFilter() {}

  virtual ~SolutionFilter() {}

  virtual SolutionFilter *clone() const = 0;

public:
  virtual bool accept(Strings const &) = 0;
};

typedef std::shared_ptr<SolutionFilter> SolutionFilterOP;

class NoExclusionFilter : public SolutionFilter {
public:
  NoExclusionFilter() : SolutionFilter() {}

  ~NoExclusionFilter() {}

  SolutionFilter *clone() const { return new NoExclusionFilter(*this); }

public:
  bool accept(Strings const &motif_names) { return true; }
};

// removes helices in the same position of the same length mostly useful for
// flex helices
class RemoveDuplicateHelices : public SolutionFilter {
public:
  RemoveDuplicateHelices() : SolutionFilter() { _seen = StringIntMap(); }

  ~RemoveDuplicateHelices() {}

  SolutionFilter *clone() const { return new RemoveDuplicateHelices(*this); }

public:
  bool accept(Strings const &motif_names) {
    _key = "";
    for (auto const &m_name : motif_names) {
      if (m_name.substr(0, 10) == "HELIX.FLEX") {
        auto spl = base::string::split(m_name, ".");
        _key += spl[0] + "." + spl[1] + "." + spl[2] + ";";
      } else {
        _key += m_name + ";";
      }
    }
    if (_seen.find(_key) != _seen.end()) {
      return false;
    }
    _seen[_key] = 1;

    return true;
  }

private:
  StringIntMap _seen;
  String _key;
};

} // namespace motif_search

#endif // RNAMAKE_NEW_SOLUTION_FILTER_H
