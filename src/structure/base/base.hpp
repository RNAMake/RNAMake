//
// Created by Joe Yesselman on 6/25/22.
//

#ifndef RNAMAKE_SRC_STRUCTURE_BASE_HPP_
#define RNAMAKE_SRC_STRUCTURE_BASE_HPP_

#include <base/string.hpp>
#include <math/rotation.hpp>
#include <util/motif_type.h>
#include <util/uuid.h>

namespace structure::base {

typedef size_t Cutpoint;
typedef std::vector<Cutpoint> Cutpoints;

class StructureException : public std::runtime_error {
public:
  /**
   * Standard constructor for StructureException
   * @param   message   Error message for structure
   */
  explicit StructureException(String const &message)
      : std::runtime_error(message) {}
};

enum class BasepairType { WC, GU, NC };

enum class ResidueType { RNA, PROTEIN, OTHER };

template <typename Restype>
String generate_bp_name(Restype const &res1, Restype const &res2) {

  auto res1_name = String("");
  auto res2_name = String("");

  if (res1.get_i_code() == ' ') {
    res1_name = res1.get_chain_id() + std::to_string(res1.get_num());
  } else {
    res1_name = res1.get_chain_id() + std::to_string(res1.get_num()) +
                res1.get_i_code();
  }

  if (res2.get_i_code() == ' ') {
    res2_name = res2.get_chain_id() + std::to_string(res2.get_num());
  } else {
    res2_name = res2.get_chain_id() + std::to_string(res2.get_num()) +
                res2.get_i_code();
  }

  if (res1.get_chain_id() < res2.get_chain_id()) {
    return res1_name + "-" + res2_name;
  }
  if (res2.get_chain_id() < res1.get_chain_id()) {
    return res2_name + "-" + res1_name;
  }

  if (res1.get_num() < res2.get_num()) {
    return res1_name + "-" + res2_name;
  } else {
    return res2_name + "-" + res1_name;
  }
}

/*template <typename Basepair, typename Structure>
base::VectorContainerOP<Index>
get_end_indexes_from_basepairs(Structure const &s,
                               std::vector<Basepair> const &bps) {
  auto start_chain_end_uuids = std::vector<util::Uuid>();
  auto end_chain_end_uuids = std::vector<util::Uuid>();

  for (auto const &r : s) {
    if (s.is_residue_start_of_chain(r)) {
      start_chain_end_uuids.push_back(r.get_uuid());
    }
    if (s.is_residue_end_of_chain(r)) {
      end_chain_end_uuids.push_back(r.get_uuid());
    }
  }

  auto end_indexes = Indexes();
  auto i = -1;
  for (auto const &bp : bps) {
    i++;
    if (bp.get_bp_type() == BasepairType::NC) {
      continue;
    }

    if (std::find(start_chain_end_uuids.begin(), start_chain_end_uuids.end(),
                  bp.get_res1_uuid()) != start_chain_end_uuids.end() &&
        std::find(end_chain_end_uuids.begin(), end_chain_end_uuids.end(),
                  bp.get_res2_uuid()) != end_chain_end_uuids.end()) {
      end_indexes.push_back(i);
    } else if (std::find(start_chain_end_uuids.begin(),
                         start_chain_end_uuids.end(),
                         bp.get_res2_uuid()) != start_chain_end_uuids.end() &&
               std::find(end_chain_end_uuids.begin(), end_chain_end_uuids.end(),
                         bp.get_res1_uuid()) != end_chain_end_uuids.end()) {
      end_indexes.push_back(i);
    }
  }

  return std::make_shared<base::VectorContainer<Index>>(end_indexes);
}     */

template <typename BPtype, typename Restype>
BPtype const *get_res_wc_or_gu_basepair(std::vector<BPtype> const &basepairs,
                                        Restype const &r) {

  for (auto const &bp : basepairs) {
    if (bp.get_bp_type() == BasepairType::NC) {
      continue;
    }
    if (bp.get_res1_uuid() == r.get_uuid() ||
        bp.get_res2_uuid() == r.get_uuid()) {
      return &bp;
    }
  }
  return nullptr;
}

template <typename Structuretype, typename Chaintype, typename BPtype,
          typename Restype>
String generate_end_id(Structuretype const &s, std::vector<BPtype> const &bps,
                       BPtype const &end) {

  auto open_chains = std::vector<Chaintype const *>();
  auto chains = s.get_chains();
  for (auto const &c : *chains) {
    if (c.get_first().get_uuid() == end.get_res1_uuid() ||
        c.get_first().get_uuid() == end.get_res2_uuid()) {
      open_chains.push_back(&c);
      break;
    }
  }

  if (open_chains.size() == 0) {
    throw std::runtime_error("could not find chain to start with");
  }

  auto seen_res = std::map<util::Uuid, int>();
  auto seen_bps = std::map<BPtype const *, int>();
  auto seen_chains = std::map<Chaintype const *, int>();
  seen_chains[open_chains[0]] = 1;

  BPtype const *bp = nullptr;
  auto ss_chains = std::vector<Strings>();
  auto seq = String("");
  auto ss = String("");
  auto dot_bracket = ' ';

  auto best_chains = std::vector<Chaintype const *>();
  Chaintype const *best_chain;
  Chaintype const *c = nullptr;
  auto best_score = 0;
  auto score = 0;
  auto pos = 0;

  while (open_chains.size() > 0) {
    c = open_chains[0];
    open_chains.erase(open_chains.begin());

    for (auto const &r : *c) {
      dot_bracket = '.';
      bp = get_res_wc_or_gu_basepair(bps, r);
      if (bp != nullptr && bp->get_bp_type() != BasepairType::NC) {
        auto &partner_res_uuid = bp->get_partner(r.get_uuid());
        auto &partner_res = s.get_residue(partner_res_uuid);
        if (seen_bps.find(bp) == seen_bps.end() &&
            seen_res.find(r.get_uuid()) == seen_res.end() &&
            seen_res.find(partner_res.get_uuid()) == seen_res.end()) {
          seen_res[r.get_uuid()] = 1;
          dot_bracket = '(';
        } else if (seen_res.find(partner_res.get_uuid()) != seen_res.end()) {
          if (seen_res[partner_res.get_uuid()] > 1) {
            dot_bracket = '.';
          } else {
            dot_bracket = ')';
            seen_res[r.get_uuid()] = 1;
            seen_res[partner_res.get_uuid()] += 1;
          }
        }
      }
      ss += dot_bracket;
      seq += r.get_name();
      if (bp != nullptr) {
        seen_bps[bp] = 1;
      }
    }

    auto dummy_str = Strings(2);
    dummy_str[0] = seq;
    dummy_str[1] = ss;
    ss_chains.push_back(dummy_str);
    ss = "";
    seq = "";
    best_score = -1;
    best_chains = std::vector<Chaintype const *>();
    for (auto const &c : *chains) {
      if (seen_chains.find(&c) != seen_chains.end()) {
        continue;
      }
      score = 0;
      for (auto const &r : c) {
        bp = get_res_wc_or_gu_basepair(bps, r);
        if (bp != nullptr && seen_bps.find(bp) == seen_bps.end()) {
          score += 1;
        }
      }
      if (score > best_score) {
        best_score = score;
      }
    }

    for (auto const &c : *chains) {
      if (seen_chains.find(&c) != seen_chains.end()) {
        continue;
      }
      score = 0;
      for (auto const &r : c) {
        bp = get_res_wc_or_gu_basepair(bps, r);
        if (bp != nullptr && seen_bps.find(bp) == seen_bps.end()) {
          score += 1;
        }
      }
      if (score == best_score) {
        best_chains.push_back(&c);
      }
    }

    best_chain = nullptr;
    best_score = 1000;
    for (auto const &c : best_chains) {
      auto i = -1;
      pos = 1000;
      for (auto const &r : *c) {
        i++;
        bp = get_res_wc_or_gu_basepair(bps, r);
        if (bp != nullptr && seen_bps.find(bp) != seen_bps.end()) {
          pos = i;
          break;
        }
      }
      if (pos < best_score) {
        best_score = pos;
        best_chain = c;
      }
    }

    if (best_chain == nullptr) {
      break;
    }
    seen_chains[best_chain] = 1;
    open_chains.push_back(best_chain);
  }

  auto ss_id = String("");
  auto i = 0;
  for (auto const &ss_chain : ss_chains) {
    ss_id += ss_chain[0] + "_";
    for (auto const &e : ss_chain[1]) {
      if (e == '(') {
        ss_id += "L";
      } else if (e == ')') {
        ss_id += "R";
      } else if (e == '.') {
        ss_id += "U";
      } else {
        throw StructureException("unexpected symbol in dot bracket notation:"
                                 " " +
                                 std::to_string(e));
      }
    }
    if (i != ss_chains.size() - 1) {
      ss_id += "_";
    }
    i++;
  }

  return ss_id;
};

template <typename Structuretype, typename Chaintype, typename BPtype,
          typename Restype>
String generate_secondary_structure(Structuretype const &s,
                                    std::vector<BPtype> const &bps) {

  auto open_chains = std::vector<Chaintype const *>();
  auto chains = s.get_chains();

  if (chains->size() == 0) {
    return String("");
  }
  open_chains.push_back(&chains->get_data()[0]);

  auto seen_res = std::map<util::Uuid, int>();
  auto seen_bps = std::map<BPtype const *, int>();
  auto seen_chains = std::map<Chaintype const *, int>();
  seen_chains[open_chains[0]] = 1;

  BPtype const *bp = nullptr;
  auto seq = String("");
  auto ss = String("");
  auto dot_bracket = ' ';

  auto best_chains = std::vector<Chaintype const *>();
  Chaintype const *best_chain;
  Chaintype const *c = nullptr;
  auto best_score = 0;
  auto score = 0;
  auto pos = 0;
  auto final_ss = String("");

  while (open_chains.size() > 0) {
    c = open_chains[0];
    open_chains.erase(open_chains.begin());

    for (auto const &r : *c) {
      dot_bracket = '.';
      bp = get_res_wc_or_gu_basepair(bps, r);
      if (bp != nullptr && bp->get_bp_type() != BasepairType::NC) {
        auto &partner_res_uuid = bp->get_partner(r.get_uuid());
        auto &partner_res = s.get_residue(partner_res_uuid);
        if (seen_bps.find(bp) == seen_bps.end() &&
            seen_res.find(r.get_uuid()) == seen_res.end() &&
            seen_res.find(partner_res.get_uuid()) == seen_res.end()) {
          seen_res[r.get_uuid()] = 1;
          dot_bracket = '(';
        } else if (seen_res.find(partner_res.get_uuid()) != seen_res.end()) {
          if (seen_res[partner_res.get_uuid()] > 1) {
            dot_bracket = '.';
          } else {
            dot_bracket = ')';
            seen_res[r.get_uuid()] = 1;
            seen_res[partner_res.get_uuid()] += 1;
          }
        }
      }
      ss += dot_bracket;
      seq += r.get_name();
      if (bp != nullptr) {
        seen_bps[bp] = 1;
      }
    }

    if (final_ss.length() > 0) {
      final_ss += "&";
    }
    final_ss += ss;

    ss = "";
    seq = "";
    best_score = -1;
    best_chains = std::vector<Chaintype const *>();
    for (auto const &c : *chains) {
      if (seen_chains.find(&c) != seen_chains.end()) {
        continue;
      }
      score = 0;
      for (auto const &r : c) {
        bp = get_res_wc_or_gu_basepair(bps, r);
        if (bp != nullptr && seen_bps.find(bp) == seen_bps.end()) {
          score += 1;
        }
      }
      if (score > best_score) {
        best_score = score;
      }
    }

    for (auto const &c : *chains) {
      if (seen_chains.find(&c) != seen_chains.end()) {
        continue;
      }
      score = 0;
      for (auto const &r : c) {
        bp = get_res_wc_or_gu_basepair(bps, r);
        if (bp != nullptr && seen_bps.find(bp) == seen_bps.end()) {
          score += 1;
        }
      }
      if (score == best_score) {
        best_chains.push_back(&c);
      }
    }

    best_chain = nullptr;
    best_score = 1000;
    for (auto const &c : best_chains) {
      auto i = -1;
      pos = 1000;
      for (auto const &r : *c) {
        i++;
        bp = get_res_wc_or_gu_basepair(bps, r);
        if (bp != nullptr && seen_bps.find(bp) != seen_bps.end()) {
          pos = i;
          break;
        }
      }
      if (pos < best_score) {
        best_score = pos;
        best_chain = c;
      }
    }

    if (best_chain == nullptr) {
      break;
    }
    seen_chains[best_chain] = 1;
    open_chains.push_back(best_chain);
  }

  return final_ss;
};

String get_dot_bracket_from_end_id(String const &);

} // namespace structure::base
#endif // RNAMAKE_SRC_STRUCTURE_BASE_HPP_
