//
// Created by Joe Yesselman on 6/30/22.
//

#ifndef RNAMAKE_SRC_STRUCTURE_BASE_POSE_HPP_
#define RNAMAKE_SRC_STRUCTURE_BASE_POSE_HPP_

#include <structure/base/structure.hpp>

namespace structure::base {

template <typename Basepair, typename Structure, typename Chain,
          typename Residue>
class Pose {
public: // types //////////////////////////////////////////////////////////////
  typedef std::vector<Residue> Residues;
  typedef std::vector<Basepair> Basepairs;

public:
  Pose(Structure &structure, Structure &proteins, Structure &small_molecules,
       Basepairs &basepairs, Indexes &end_indexes, Strings &end_ids,
       String &name, String &dot_bracket)
      : _structure(std::move(structure)), _proteins(std::move(proteins)),
        _small_molecules(std::move(small_molecules)),
        _basepairs(std::move(basepairs)), _end_indexes(end_indexes),
        _end_ids(end_ids), _name(name), _dot_bracket(std::move(dot_bracket)) {

    /*expects<StructureException>(
        _end_ids.size() == _end_indexes.size(),
        "Pose must have the same number of ends as end_ids has " +
            std::to_string(_end_indexes.size()) + " ends and " +
            std::to_string(end_ids.size()) + "end ids");   */
  }

public: // iterators //////////////////////////////////////////////////////////
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

  inline const Cutpoints &get_cutpoints() const {
    return _structure.get_cutpoints();
  };

  inline void get_chains() const {
    // return _structure.get_chains();
  }

  inline bool is_residue_start_of_chain(const Residue &r) const {
    return _structure.is_residue_start_of_chain(r);
  }

  inline bool is_residue_end_of_chain(const Residue &r) const {
    return _structure.is_residue_end_of_chain(r);
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

    /* std::cout << index << " " << _basepairs.size() << " " << _end_indexes
                                                                    .size()
               << std::endl;
     std::cout << _end_indexes[index] << std::endl;
     std::cout << _basepairs[_end_indexes[index]].get_name_str() << std::endl;
      */
    return _basepairs[_end_indexes[index]];
  }

public: // only available to ends with ref frames?
  inline const math::Matrix3x3 &get_end_ref_frame(Index end_index) const {
    return get_end(end_index).get_ref_frame();
  }

  inline const math::Vector3 &get_end_center(Index end_index) const {
    return get_end(end_index).get_center();
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

public: // trivial getters ////////////////////////////////////////////////////
  inline const String &get_dot_bracket() const { return _dot_bracket; }

  inline const Indexes &get_end_indexes() const { return _end_indexes; }

  inline const Strings &get_end_ids() const { return _end_ids; }

  inline const String &get_name() const { return _name; }

public: // other getters
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

public:
  void move(const math::Vector3 &p) {
    _structure.move(p);
    for (auto &bp : _basepairs) {
      bp.move(p);
    }
  }

  void rotate(const math::Matrix3x3 rot) {
    _structure.rotate(rot);
    for (auto &bp : _basepairs) {
      bp.rotate(rot);
    }
  }

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

}


#endif // RNAMAKE_SRC_STRUCTURE_BASE_POSE_HPP_