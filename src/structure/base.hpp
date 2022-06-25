//
// Created by Joe Yesselman on 6/25/22.
//

#ifndef RNAMAKE_SRC_STRUCTURE_BASE_HPP_
#define RNAMAKE_SRC_STRUCTURE_BASE_HPP_

#include <base/string.hpp>
#include <util/uuid.h>
#include <util/motif_type.h>

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

class Basepair {
public:
  inline Basepair(const util::Uuid &res1_uuid, const util::Uuid &res2_uuid,
                  const util::Uuid &uuid, const BasepairType &bp_type,
                  String &name)
      : _res1_uuid(res1_uuid), _res2_uuid(res2_uuid), _uuid(uuid),
        _bp_type(bp_type), _name(std::move(name)) {}

  inline Basepair(Basepair const &bp) = default;

  virtual ~Basepair() = default;

public:
  inline bool operator==(Basepair const &other) const {
    return _uuid == other._uuid;
  }

  inline bool operator!=(Basepair const &other) const {
    return _uuid != other._uuid;
  }

public:
  [[nodiscard]] util::Uuid const &get_partner(util::Uuid const &uuid) const {
    if (uuid == _res1_uuid) {
      return _res2_uuid;
    } else {
      return _res1_uuid;
    }
  }

  [[nodiscard]] inline BasepairType const &get_bp_type() const {
    return _bp_type;
  }

  [[nodiscard]] inline util::Uuid const &get_uuid() const { return _uuid; }

  [[nodiscard]] inline const String &get_name() const { return _name; }

  [[nodiscard]] inline util::Uuid const &get_res1_uuid() const {
    return _res1_uuid;
  }

  [[nodiscard]] inline util::Uuid const &get_res2_uuid() const {
    return _res2_uuid;
  }

private:
  util::Uuid _uuid;
  util::Uuid _res1_uuid;
  util::Uuid _res2_uuid;
  BasepairType _bp_type;
  String _name;
};

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
  inline Structure(Residues &res, Cutpoints &cut_points)
      : _residues(std::move(res)), _cut_points(std::move(cut_points)) {}

  ~Structure() = default;

public: // res iterator
  typedef typename Residues::const_iterator const_iterator;

  const_iterator begin() const { return _residues.begin(); }
  const_iterator end() const { return _residues.end(); }

public: // get_residue interface
  Restype const &get_residue(int num, char chain_id, char i_code) const {

    for (auto const &r : _residues) {
      if (num == r.get_num() && chain_id == r.get_chain_id() &&
          i_code == r.get_i_code()) {
        return r;
      }
    }

    /*auto ss = std::stringstream();
    ss << "cannot find residue with num: " << num << " chain id: " << chain_id
       << " and i_code";
    throw StructureException(ss.str());      */
  }

  // TODO Need to remove this function and only use char
  Restype const &get_residue(int num, String chain_id, String i_code) const {

    for (auto const &r : _residues) {
      if (num == r.get_num() && chain_id[0] == r.get_chain_id() &&
          i_code[0] == r.get_i_code()) {
        return r;
      }
    }

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
  Pose(Structure const &structure, std::vector<Basepair> const &basepairs,
       Indexes const &end_indexes, Strings const &end_ids, String const &name)
      : structure_(structure), basepairs_(basepairs), end_indexes_(end_indexes),
        end_ids_(end_ids), name_(name) {

    /*expects<StructureException>(
        end_ids_.size() == end_indexes_.size(),
        "Pose must have the same number of ends as end_ids has " +
            std::to_string(end_indexes_.size()) + " ends and " +
            std::to_string(end_ids.size()) + "end ids");   */
  }

protected:
  // let dervived classes fill in members
  Pose() : structure_(Structure(Residues(), Cutpoints())) {}

public: // iterators
  // residue iterator
  typedef typename Residues::const_iterator const_iterator;

  const_iterator begin() const { return structure_.begin(); }
  const_iterator end() const { return structure_.end(); }

  // basepair iterator
  typedef typename std::vector<Basepair>::const_iterator const_bp_iterator;

  const_bp_iterator bps_begin() const { return basepairs_.begin(); }
  const_bp_iterator bps_end() const { return basepairs_.end(); }

  // end iterator
  Indexes::const_iterator end_indexes_begin() const {
    return end_indexes_.begin();
  }
  Indexes::const_iterator end_indexes_end() const { return end_indexes_.end(); }

public: // structure wrappers
  inline String get_sequence() { return structure_.get_sequence(); }

  inline Residue const &get_residue(int num, char chain_id, char i_code) const {
    return structure_.get_residue(num, chain_id, i_code);
  }

  inline Residue const &get_residue(util::Uuid const &uuid) const {
    return structure_.get_residue(uuid);
  }

  inline Residue const &get_residue(Index index) const {
    return structure_.get_residue(index);
  }

  inline size_t get_num_residues() const {
    return structure_.get_num_residues();
  }

  inline size_t get_num_chains() const { return structure_.get_num_chains(); }

  inline void get_chains() const {
    // return structure_.get_chains();
  }

public: // get basepairs interface
  /*BasepairsOP get_basepairs(util::Uuid const &bp_uuid) const {
    auto bps = std::vector<Basepair>();
    for (auto const &bp : basepairs_) {
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
    for (auto const &bp : basepairs_) {
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
    for (auto const &bp : basepairs_) {
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
    for (auto const &bp : basepairs_) {
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
    for (auto const &bp : basepairs_) {
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
    for (auto const &bp : basepairs_) {
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
    /*expects<StructureException>(index < basepairs_.size(),
                           "cannot get basepair " + std::to_string(index) +
                               " only " + std::to_string(basepairs_.size()) +
                               " total residues");   */
    return basepairs_[index];
  }

public: // get end interace
  Basepair const &get_end(util::Uuid const &bp_uuid) const {
    auto bps = std::vector<Basepair const *>();
    for (auto const &ei : end_indexes_) {
      auto &bp = basepairs_[ei];
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
    for (auto const &ei : end_indexes_) {
      auto &bp = basepairs_[ei];
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
    for (auto const &ei : end_indexes_) {
      auto &bp = basepairs_[ei];
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

    /*expects<StructureException>(index < end_indexes_.size(),
                                "trying to get end: " + std::to_string(index) +
                                    " there are only " +
                                    std::to_string(end_indexes_.size()));    */

    std::cout << index << " " << basepairs_.size() << " " << end_indexes_.size()
              << std::endl;
    std::cout << end_indexes_[index] << std::endl;
    std::cout << basepairs_[end_indexes_[index]].get_name_str() << std::endl;

    return basepairs_[end_indexes_[index]];
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
    for (auto const &end_id : end_ids_) {
      i++;
      if (end_id == nend_id) {
        bps.push_back(&basepairs_[end_indexes_[i]]);
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
    if (index >= end_ids_.size()) {
      throw StructureException(
          "trying to get end_id: " + std::to_string(index) +
          " there are only " + std::to_string(end_ids_.size()));
    }
    return end_ids_[index];
  }

  int get_end_index(base::SimpleStringCOP name) const {
    auto &bp = get_end(name);
    int i = 0;
    for (auto const &ei : end_indexes_) {
      auto &end = basepairs_[ei];
      if (bp == end) {
        return i;
      }
      i++;
    }
    throw StructureException("cannot find end with name: " + name->get_str());
  } */

  int get_end_index(String const &str) const {
    int i = 0;
    for (auto const &ei : end_ids_) {
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

  size_t get_num_basepairs() const { return basepairs_.size(); }

  size_t get_num_ends() const { return end_indexes_.size(); }

  //base::SimpleStringCOP get_name() { return name_; }

private:
  Structure structure_;
  std::vector<Basepair> basepairs_;
  Indexes end_indexes_;
  String name_;
  mutable Strings end_ids_;
};

template <typename BPtype, typename Structuretype, typename Chaintype,
          typename Restype>
class Segment : public Pose<BPtype, Structuretype, Chaintype, Restype> {
private:
  typedef Pose<BPtype, Structuretype, Chaintype, Restype> BaseClass;
  typedef std::vector<Restype> Residues;

public:
  inline Segment(Structuretype const &structure,
                 std::vector<BPtype> const &basepairs,
                 Indexes const &end_indexes,
                 String const &end_ids,
                 String  name, util::MotifType segment_type,
                 Index aligned_end_index, util::Uuid const &uuid)
      : BaseClass(structure, basepairs, end_indexes, end_ids, name),
        segment_type_(segment_type), aligned_end_index_(aligned_end_index) {}

  // TODO Change back to private
public:
  // let dervived classes fill in members
  Segment() : BaseClass() {}

  typedef std::vector<BPtype> Basepairs;

public:
  inline Index get_aligned_end_indx() const { return aligned_end_index_; }

  util::Uuid const &get_uuid() const { return uuid_; }

  util::MotifType get_segment_type() const { return segment_type_; }

protected:
  util::MotifType segment_type_;
  Index aligned_end_index_;
  util::Uuid uuid_;
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

/*template <typename BPtype, typename Structuretype>
base::VectorContainerOP<Index>
get_end_indexes_from_basepairs(Structuretype const &s,
                               std::vector<BPtype> const &bps) {
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
