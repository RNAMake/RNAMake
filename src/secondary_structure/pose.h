//
// Created by Joseph Yesselman on 10/27/17.
//

#ifndef RNAMAKE_NEW_RNA_STRUCTURE_H
#define RNAMAKE_NEW_RNA_STRUCTURE_H

#include <primitives/pose.h>
#include <secondary_structure/basepair.h>
#include <secondary_structure/chain.h>
#include <secondary_structure/residue.h>
#include <secondary_structure/structure.h>

namespace secondary_structure {

class Pose : public primitives::Pose<Basepair, Structure, Chain, Residue> {
public:
  typedef primitives::Pose<Basepair, Structure, Chain, Residue> BaseClass;

public:
  inline Pose(Structure const &s, std::vector<Basepair> const &basepairs,
              Indexes const &end_indexes, base::SimpleStringCOPs const &end_ids,
              base::SimpleStringCOP name)
      : BaseClass(s, basepairs, end_indexes, end_ids, name) {}

  inline Pose(Pose const &rs) : BaseClass(rs) {}

public:
  bool operator==(Pose const &rs) const {

    if (structure_ != rs.structure_) {
      return false;
    }
    if (basepairs_.size() != rs.basepairs_.size()) {
      return false;
    }
    if (*name_ != *rs.name_) {
      return false;
    }
    if (end_indexes_.size() != rs.end_indexes_.size()) {
      return false;
    }
    if (end_ids_.size() != rs.end_ids_.size()) {
      return false;
    }

    for (auto i = 0; i < basepairs_.size(); i++) {
      if (basepairs_[i] != rs.basepairs_[i]) {
        return false;
      }
    }
    for (auto i = 0; i < end_ids_.size(); i++) {
      if (*end_ids_[i] != *rs.end_ids_[i]) {
        return false;
      }
    }
    return true;
  }

  inline bool operator!=(Pose const &rs) const { return !(*this == rs); }

public:
  inline String get_dot_bracket() const { return structure_.get_dot_bracket(); }

public: // setters
  void set_sequence(String const &sequence) {
    structure_.set_sequence(sequence);
  }

  inline void set_residue_identity(Index residue_index, char name) {
    structure_.set_residue_identity(residue_index, name);
  }

private:
  /*std::unique_ptr<String>
  _get_basepair_str(
          Basepair const & bp) {
      auto s = std::unique_ptr<String>(new String(""));
      auto res = get_bp_res(bp);
      *s += bp.get_str() + ";";
      *s += std::to_string(res->at(0).get_num()) + "|" +
  res->at(0).get_chain_id() + "|" + res->at(0).get_i_code() +
            ";";
      *s += std::to_string(res->at(1).get_num()) + "|" +
  res->at(1).get_chain_id() + "|" + res->at(1).get_i_code() +
            "@";
      return s;
  }
  Basepair *
  _bp_from_str(
          String const & s) {
      auto bp_spl = base::split_str_by_delimiter(s, ";");
      expects<PoseException>(
              bp_spl.size() == 3,
              "incorrect number of sections in basepair string");
      auto r1_info = base::split_str_by_delimiter(bp_spl[1], "|");
      auto r2_info = base::split_str_by_delimiter(bp_spl[2], "|");
      auto & res1 = structure_.get_residue(std::stoi(r1_info[0]), r1_info[1][0],
  r1_info[2][0]); auto & res2 = structure_.get_residue(std::stoi(r2_info[0]),
  r2_info[1][0], r2_info[2][0]); return new Basepair(bp_spl[0], util::Uuid(),
  res1.get_uuid(), res2.get_uuid());
  }*/
};

typedef std::shared_ptr<Pose> PoseOP;
typedef std::unique_ptr<Pose> PoseUP;

inline base::VectorContainerOP<Index>
get_ends_from_basepairs(Structure const &s, Basepairs const &bps) {
  return primitives::get_end_indexes_from_basepairs<Basepair, Structure>(s,
                                                                         bps);
}

inline String generate_end_id(Structure const &s, Basepairs const &bps,
                              Basepair const &end) {
  return primitives::generate_end_id<Structure, Chain, Basepair, Residue>(
      s, bps, end);
}

} // namespace secondary_structure

#endif // RNAMAKE_NEW_RNA_STRUCTURE_H