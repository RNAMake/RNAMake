//
// Created by Joseph Yesselman on 12/15/17.
//

#ifndef RNAMAKE_NEW_STRUCTURE_RNA_STRUCTURE_H
#define RNAMAKE_NEW_STRUCTURE_RNA_STRUCTURE_H


#include <primitives/pose.h>
#include <structure/structure.h>
#include <structure/basepair.h>
#include <structure/pdb_parser.h>

namespace structure {

class Pose : public primitives::Pose<Basepair, Structure, Chain, Residue> {
public:
    friend class SegmentFactory;

public:
    typedef primitives::Pose<Basepair, Structure, Chain, Residue> BaseClass;

public:
    inline
    Pose(
            Structure const & structure,
            Basepairs const & basepairs,
            Indexes const & end_indexes,
            base::SimpleStringCOPs const & end_ids,
            base::SimpleStringCOP name,
            Structure const & proteins,
            Structure const & small_molecules,
            base::SimpleStringCOP dot_bracket):
            BaseClass(structure, basepairs, end_indexes, end_ids, name),
            proteins_(proteins),
            small_molecules_(small_molecules),
            dot_bracket_(dot_bracket) {}
    inline
    Pose(
            Pose const & p):
            BaseClass(p.structure_, p.basepairs_, p.end_indexes_, p.end_ids_, p.name_),
            proteins_(p.proteins_),
            small_molecules_(p.small_molecules_),
            dot_bracket_(p.dot_bracket_) {}


public:
    inline
    bool
    operator ==(Pose const & p) const {
        return is_equal(p);
    }

    inline
    bool
    operator !=(Pose const & p) const  {
        return !is_equal(p);
    }

public:

    const_iterator protein_begin() const { return proteins_.begin(); }
    const_iterator protein_end()   const { return proteins_.end(); }


public:
    bool
    is_equal(
            Pose const & p,
            bool check_uuid = true) const {
        if(basepairs_.size() != p.basepairs_.size()) { return false; }
        if(end_indexes_.size() != p.end_indexes_.size()) { return false; }
        if(*name_ != *p.name_) { return false; }
        if(*dot_bracket_ != *p.dot_bracket_) { return false; }
        if(! structure_.is_equal(p.structure_, check_uuid)) { return false; }
        if(! proteins_.is_equal(p.proteins_, check_uuid)) { return false; }
        if(! small_molecules_.is_equal(p.small_molecules_, check_uuid)) { return false; }
        return true;
    }

public:
public: // non const methods
    void
    move(
            math::Point const & p) {
        structure_.move(p);
        proteins_.move(p);
        small_molecules_.move(p);
        for(auto & bp : basepairs_) { bp.move(p); }
    }

    void
    transform(
            math::Matrix const & r,
            math::Vector const & t,
            math::Point & dummy) {
        structure_.transform(r, t, dummy);
        proteins_.transform(r, t, dummy);
        small_molecules_.transform(r, t, dummy);
        for(auto & bp : basepairs_) { bp.transform(r, t, dummy); }
    }

    inline
    void
    transform(
            math::Matrix const & r,
            math::Vector const & t) {
        auto dummy = math::Point();
        transform(r, t, dummy);
    }

public:
    inline
    String
    get_dot_bracket() { return dot_bracket_->get_str(); }


protected:
    Structure proteins_;
    Structure small_molecules_;
    base::SimpleStringCOP dot_bracket_;
};

typedef std::shared_ptr<Pose> PoseOP;

inline
base::VectorContainerOP<Index>
get_end_indexes_from_basepairs(
        Structure const & s,
        Basepairs const & bps) {
    return primitives::get_end_indexes_from_basepairs<Basepair, Structure>(s, bps);
}

inline
String
generate_end_id(
        Structure const & s,
        Basepairs const & bps,
        Basepair const & end) {
    return primitives::generate_end_id<Structure, Chain, Basepair, Residue>(s, bps, end);
}

inline
String
generate_secondary_structure(
        Structure const & s,
        Basepairs const & bps) {
    return primitives::generate_secondary_structure<Structure, Chain, Basepair, Residue>(s, bps);
}

PoseOP
get_pose_from_pdb(
        String const &,
        PDBParser &);

}
#endif //RNAMAKE_NEW_STRUCTURE_RNA_STRUCTURE_H
