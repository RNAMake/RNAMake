//
// Created by Joe Yesselman on 6/25/22.
//

#ifndef RNAMAKE_SRC_STRUCTURE_BASE_HPP_
#define RNAMAKE_SRC_STRUCTURE_BASE_HPP_

#include <base/string.hpp>
#include <util/motif_type.h>
#include <util/uuid.h>

namespace structure {

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

template <typename Residue> class Chain {
public:
  typedef std::vector<Residue> Residues;

public:
  inline explicit Chain(Residues const &residues) : _residues(residues) {
    if (_residues.empty()) {
      throw StructureException("cannot initiate an empty chain");
    }
  }

  ~Chain() = default;

public: // iterator ///////////////////////////////////////////////////////////
  typedef typename Residues::const_iterator const_iterator;

  const_iterator begin() const noexcept { return _residues.begin(); }
  const_iterator end() const noexcept { return _residues.end(); }

public:
  [[nodiscard]] inline size_t get_length() const {
    return (int)_residues.size();
  }

  [[nodiscard]] inline const Residue &get_first() const { return _residues[0]; }

  [[nodiscard]] inline const Residue &get_last() const {
    return _residues.back();
  }

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

private:
  Residues _residues;
};

template <typename Chaintype, typename Restype> class Structure {
public:
  typedef std::vector<Restype> Residues;
  typedef std::vector<Chaintype> Chains;
  typedef std::vector<Chaintype> ChainsOP;

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
  Restype const &get_residue(int num, const String &chain_id,
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

  Restype const &get_residue(util::Uuid const &uuid) const {

    for (auto const &r : _residues) {
      if (r.get_uuid() == uuid) {
        return r;
      }
    }

    throw StructureException("cannot find residue by uuid");
  }

  Restype const &get_residue(Index index) const {

    /*expects<StructureException>(
        index < _residues.size(),
        "cannot get residue " + std::to_string(index) + " only " +
            std::to_string(_residues.size()) + " total residues");   */

    return _residues[index];
  }

  int get_res_index(Restype const &res) const {
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
    auto chains = std::vector<Chaintype>();
    auto i = 0;
    for (auto const &r : _residues) {
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
    if (res.size() > 0) {
      chains.push_back(Chaintype(res));
    }
    return std::make_shared<Chains>(chains);
  }

  Cutpoints const &get_cutpoints() const { return _cut_points; }

  size_t get_num_residues() const { return _residues.size(); }

  size_t get_num_chains() const { return _cut_points.size(); }

  String get_sequence() const {
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

  bool is_residue_start_of_chain(Restype const &r) const {
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

  bool is_residue_end_of_chain(Restype const &r) const {
    auto res_index = get_res_index(r);
    for (auto const c : _cut_points) {
      if (res_index == c - 1) {
        return true;
      }
    }
    return false;
  }

private:
  Residues _residues;
  Cutpoints _cut_points;
};

template <typename Basepair, typename Structure, typename Chain,
          typename Residue>
class Pose {
public: // types
  typedef std::vector<Residue> Residues;
  typedef std::vector<Basepair> Basepairs;

public:
  Pose(Structure &structure, Structure &proteins, Structure &small_molecules,
       Basepairs &basepairs, Indexes &end_indexes, Strings &end_ids,
       String &name)
      : _structure(std::move(structure)), _proteins(std::move(proteins)),
        _small_molecules(std::move(small_molecules)),
        _basepairs(std::move(basepairs)), _end_indexes(end_indexes),
        _end_ids(end_ids), _name(name) {

    /*expects<StructureException>(
        _end_ids.size() == _end_indexes.size(),
        "Pose must have the same number of ends as end_ids has " +
            std::to_string(_end_indexes.size()) + " ends and " +
            std::to_string(end_ids.size()) + "end ids");   */
  }

public: // iterators
  // residue iterator
  typedef typename Residues::const_iterator const_iterator;

  const_iterator begin() const { return _structure.begin(); }
  const_iterator end() const { return _structure.end(); }

  // basepair iterator
  typedef typename std::vector<Basepair>::const_iterator const_bp_iterator;

  const_bp_iterator bps_begin() const { return _basepairs.begin(); }
  const_bp_iterator bps_end() const { return _basepairs.end(); }

  // end iterator
  Indexes::const_iterator end_indexes_begin() const {
    return _end_indexes.begin();
  }
  Indexes::const_iterator end_indexes_end() const { return _end_indexes.end(); }

public: // structure wrappers
  inline String get_sequence() { return _structure.get_sequence(); }

  inline Residue const &get_residue(int num, char chain_id, char i_code) const {
    return _structure.get_residue(num, chain_id, i_code);
  }

  inline Residue const &get_residue(util::Uuid const &uuid) const {
    return _structure.get_residue(uuid);
  }

  inline Residue const &get_residue(Index index) const {
    return _structure.get_residue(index);
  }

  inline size_t get_num_residues() const {
    return _structure.get_num_residues();
  }

  inline size_t get_num_chains() const { return _structure.get_num_chains(); }

  inline void get_chains() const {
    // return _structure.get_chains();
  }

public: // get basepairs interface
  /*BasepairsOP get_basepairs(util::Uuid const &bp_uuid) const {
    auto bps = std::vector<Basepair>();
    for (auto const &bp : _basepairs) {
      if (bp.get_uuid() == bp_uuid) {
        bps.push_back(bp);
      }
      if (bp.get_res1_uuid() == bp_uuid || bp.get_res2_uuid() == bp_uuid) {
        bps.push_back(bp);
      }
    }

    if (bps.size() == 0) {
      throw StructureException(
          "could not find any basepairs with this uuid for "
          "residues or basepairs");
    }

    return std::make_shared<base::VectorContainer<Basepair>>(bps);
  }

  BasepairsOP get_basepairs(util::Uuid const &uuid1,
                            util::Uuid const &uuid2) const {

    auto bps = std::vector<Basepair>();
    for (auto const &bp : _basepairs) {
      if (bp.get_res1_uuid() == uuid1 && bp.get_res2_uuid() == uuid2) {
        bps.push_back(bp);
      }
      if (bp.get_res1_uuid() == uuid2 && bp.get_res2_uuid() == uuid1) {
        bps.push_back(bp);
      }
    }

    if (bps.size() == 0) {
      throw StructureException("could not find any basepairs with these two
  uuids");
    }

    return std::make_shared<base::VectorContainer<Basepair>>(bps);
  }

  BasepairsOP get_basepairs(String const &name) const {
    auto bps = std::vector<Basepair>();
    for (auto const &bp : _basepairs) {
      if (name == bp.get_name()->get_str()) {
        bps.push_back(bp);
      }
    }

    if (bps.size() == 0) {
      throw StructureException("could not find any basepairs with this name: " +
                          name);
    }

    return std::make_shared<base::VectorContainer<Basepair>>(bps);
  }  */

public: // get basepair interface  (single basepair!)
  Basepair const &get_basepair(util::Uuid const &bp_uuid) const {
    auto bps = std::vector<Basepair const *>();
    for (auto const &bp : _basepairs) {
      if (bp.get_uuid() == bp_uuid) {
        bps.push_back(&bp);
      }
      if (bp.get_res1_uuid() == bp_uuid || bp.get_res2_uuid() == bp_uuid) {
        bps.push_back(&bp);
      }
    }
    if (bps.size() > 1) {
      throw StructureException("got more than one basepair matching this uuid");
    }
    if (bps.size() == 1) {
      return *bps[0];
    } else {
      throw StructureException("no basepairs matching this uuid");
    }
  }

  Basepair const &get_basepair(util::Uuid const &uuid1,
                               util::Uuid const &uuid2) const {

    auto bps = std::vector<Basepair const *>();
    for (auto const &bp : _basepairs) {
      if (bp.get_res1_uuid() == uuid1 && bp.get_res2_uuid() == uuid2) {
        bps.push_back(&bp);
      }
      if (bp.get_res1_uuid() == uuid2 && bp.get_res2_uuid() == uuid1) {
        bps.push_back(&bp);
      }
    }
    if (bps.size() > 1) {
      throw StructureException(
          "got more than one basepair matching residue uuids");
    }
    if (bps.size() == 1) {
      return *bps[0];
    } else {
      throw StructureException(
          "no basepair found matching residue uuids supplied");
    }
  }

  Basepair const &get_basepair(String const &name) const {
    auto bps = std::vector<Basepair const *>();
    for (auto const &bp : _basepairs) {
      if (bp.get_name()->get_str() == name) {
        bps.push_back(&bp);
      }
    }

    if (bps.size() > 1) {
      throw StructureException(
          "got more than one basepair matching this name: " + name);
    }
    if (bps.size() == 1) {
      return *bps[0];
    } else {
      throw StructureException(
          "no basepair found matching residue uuids supplied");
    }
  }

  inline Basepair const &get_basepair(Index index) const {
    /*expects<StructureException>(index < _basepairs.size(),
                           "cannot get basepair " + std::to_string(index) +
                               " only " + std::to_string(_basepairs.size()) +
                               " total residues");   */
    return _basepairs[index];
  }

public: // get end interace
  Basepair const &get_end(util::Uuid const &bp_uuid) const {
    auto bps = std::vector<Basepair const *>();
    for (auto const &ei : _end_indexes) {
      auto &bp = _basepairs[ei];
      if (bp.get_uuid() == bp_uuid) {
        bps.push_back(&bp);
      }
      if (bp.get_res1_uuid() == bp_uuid || bp.get_res2_uuid() == bp_uuid) {
        bps.push_back(&bp);
      }
    }
    if (bps.size() > 1) {
      throw StructureException("got more than one basepair matching this uuid");
    }
    if (bps.size() == 1) {
      return *bps[0];
    } else {
      throw StructureException("no end found matching basepair uuid supplied");
    }
  }

  Basepair const &get_end(util::Uuid const &uuid1,
                          util::Uuid const &uuid2) const {

    auto bps = std::vector<Basepair const *>();
    for (auto const &ei : _end_indexes) {
      auto &bp = _basepairs[ei];
      if (bp.get_res1_uuid() == uuid1 && bp.get_res2_uuid() == uuid2) {
        bps.push_back(&bp);
      }
      if (bp.get_res1_uuid() == uuid2 && bp.get_res2_uuid() == uuid1) {
        bps.push_back(&bp);
      }
    }
    if (bps.size() > 1) {
      throw StructureException("got more than one end matching residue uuids");
    }
    if (bps.size() == 1) {
      return *bps[0];
    } else {
      throw StructureException("no end found matching residue uuids supplied");
    }
  }

  Basepair const &get_end(String const &name) const {
    auto bps = std::vector<Basepair const *>();
    for (auto const &ei : _end_indexes) {
      auto &bp = _basepairs[ei];
      if (bp.get_name()->get_str() == name) {
        bps.push_back(&bp);
      }
    }

    if (bps.size() > 1) {
      throw StructureException(
          "got more than one basepair matching this name: " + name);
    }
    if (bps.size() == 1) {
      return *bps[0];
    } else {
      throw StructureException("cannot find end with name: " + name);
    }
  }

  inline Basepair const &get_end(Index index) const {

    /*expects<StructureException>(index < _end_indexes.size(),
                                "trying to get end: " + std::to_string(index) +
                                    " there are only " +
                                    std::to_string(_end_indexes.size()));    */

    std::cout << index << " " << _basepairs.size() << " " << _end_indexes.size()
              << std::endl;
    std::cout << _end_indexes[index] << std::endl;
    std::cout << _basepairs[_end_indexes[index]].get_name_str() << std::endl;

    return _basepairs[_end_indexes[index]];
  }

public:
  inline String get_end_name(Index index) const {
    return get_end(index).get_name_str();
  }

public: // get end by end id
  // avoid confliction with getting by name ... not pretty
  Basepair const &get_end_by_id(String const &nend_id) const {
    auto bps = std::vector<Basepair const *>();
    int i = -1;
    for (auto const &end_id : _end_ids) {
      i++;
      if (end_id == nend_id) {
        bps.push_back(&_basepairs[_end_indexes[i]]);
      }
    }

    if (bps.size() > 1) {
      throw StructureException(
          "got more than one basepair matching this end_id: " + nend_id);
    }
    if (bps.size() == 1) {
      return *bps[0];
    } else {
      throw StructureException("cannot find end with end_id: " + nend_id);
    }
  }

public: // other getters
  /*base::SimpleStringCOP get_end_id(Index index) const {
    if (index >= _end_ids.size()) {
      throw StructureException(
          "trying to get end_id: " + std::to_string(index) +
          " there are only " + std::to_string(_end_ids.size()));
    }
    return _end_ids[index];
  }

  int get_end_index(base::SimpleStringCOP name) const {
    auto &bp = get_end(name);
    int i = 0;
    for (auto const &ei : _end_indexes) {
      auto &end = _basepairs[ei];
      if (bp == end) {
        return i;
      }
      i++;
    }
    throw StructureException("cannot find end with name: " + name->get_str());
  } */

  int get_end_index(String const &str) const {
    int i = 0;
    for (auto const &ei : _end_ids) {
      if (ei == str) {
        return i;
      }
      if (get_end(i).get_name()->get_str() == str) {
        return i;
      }
      i++;
    }
    throw StructureException("cannot find end with str: " + str);
  }

  /*ResiduesOP get_bp_res(Basepair const &bp) const {
    auto res = std::vector<Residue>();
    res.push_back(get_residue(bp.get_res1_uuid()));
    res.push_back(get_residue(bp.get_res2_uuid()));
    return std::make_shared<base::VectorContainer<Residue>>(res);
  } */

  size_t get_num_basepairs() const { return _basepairs.size(); }

  size_t get_num_ends() const { return _end_indexes.size(); }

  // base::SimpleStringCOP get_name() { return _name; }

protected:
  Structure _structure;
  Structure _proteins;
  Structure _small_molecules;
  Basepairs _basepairs;
  Indexes _end_indexes;
  String _name;
  String _dot_bracket;
  mutable Strings _end_ids;
};

template <typename Basepair, typename Structure, typename Chain,
          typename Residue>
class Segment : public Pose<Basepair, Structure, Chain, Residue> {
private:
  typedef Pose<Basepair, Structure, Chain, Residue> BaseClass;
  typedef std::vector<Residue> Residues;
  typedef std::vector<Basepair> Basepairs;

public:
  inline Segment(Structure &structure, Structure &proteins,
                 Structure &small_molecules, Basepairs &basepairs,
                 Indexes &end_indexes, Strings &end_ids, String &name,
                 util::MotifType segment_type, Index aligned_end_index,
                 util::Uuid const &uuid)
      : BaseClass(structure, proteins, small_molecules, basepairs, end_indexes,
                  end_ids, name),
        _segment_type(segment_type), _aligned_end_index(aligned_end_index),
        _uuid(uuid) {}

public:
  [[nodiscard]] inline Index get_aligned_end_indx() const {
    return _aligned_end_index;
  }

  [[nodiscard]] util::Uuid const &get_uuid() const { return _uuid; }

  [[nodiscard]] util::MotifType get_segment_type() const {
    return _segment_type;
  }

private:
  util::MotifType _segment_type;
  Index _aligned_end_index;
  util::Uuid _uuid;
};

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

} // namespace structure
#endif // RNAMAKE_SRC_STRUCTURE_BASE_HPP_
