#ifndef PYBIND11_MOTIF_HPP
#define PYBIND11_MOTIF_HPP

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/cast.h>
#include <pybind11/operators.h>
#include <memory>

// motif includes
#include <motif/motif.h>
#include <motif/motif_ensemble.h>
#include <motif/motif_factory.h>
#include <motif/motif_scorer.h>
#include <motif/motif_state.h>
#include <motif/motif_state_aligner.h>
#include <motif/motif_state_ensemble.h>
#include <motif/motif_to_secondary_structure.h>
#include <motif/pose.h>
#include <motif/pose_factory.h>

namespace motif {
    namespace py = pybind11;

    void
    add_bindings(py::module_ & m ) {
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// motif
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // free functions

        m.def("align_motif", [] (structure::BasepairStateOP const & ref_bp_state, structure::BasepairOP const & motif_end, motif::MotifOP & motif) {
                  align_motif(ref_bp_state, motif_end, motif); },
              py::arg("ref_bp_state"),
              py::arg("motif_end"),
              py::arg("motif")
        );

        m.def("align_motif_state", [] (structure::BasepairStateOP const & ref_bp_state, motif::MotifStateOP org_state) {
                  align_motif_state(ref_bp_state, org_state); },
              py::arg("ref_bp_state"),
              py::arg("org_state")
        );

        m.def("clash_between_motifs", [] (motif::MotifOP m1, motif::MotifOP m2, double clash_radius) ->  int {
                  return clash_between_motifs(m1, m2, clash_radius); },
              py::arg("m1"),
              py::arg("m2"),
              py::arg("clash_radius") = 2.7
        );

        m.def("file_to_motif", [] (String const & path) -> motif::MotifOP {
                  return file_to_motif(path); },
              py::arg("path")
        );

        m.def("get_aligned_motif", [] (structure::BasepairOP const & ref_bp, structure::BasepairOP const & motif_end, motif::MotifOP const & motif) -> motif::MotifOP {
                  return get_aligned_motif(ref_bp, motif_end, motif); },
              py::arg("ref_bp"),
              py::arg("motif_end"),
              py::arg("motif")
        );

        m.def("get_aligned_motif_state", [] (structure::BasepairStateOP const & ref_bp_state, motif::MotifStateOP cur_state, motif::MotifStateOP org_state) {
                  get_aligned_motif_state(ref_bp_state, cur_state, org_state); },
              py::arg("ref_bp_state"),
              py::arg("cur_state"),
              py::arg("org_state")
        );

        m.def("ref_motif", [] () -> motif::Motif {
            return ref_motif(); }
        );
        // classes

        py::class_<motif::Motif, std::shared_ptr<motif::Motif>>(m, "Motif")
                // ctors
                .def(py::init<>())
                .def(py::init<structure::StructureOP const &,structure::BasepairOPs const &,structure::BasepairOPs const &>())
                .def(py::init<String const,structure::ResidueTypeSet const>())
                .def(py::init<motif::Motif const &>())
                .def(py::init<structure::RNAStructure const &>())
                        // methods
                .def("transform",[] (motif::Motif  & ptr, math::Transform const & t) {
                    ptr.transform(t); } )
                .def("move",[] (motif::Motif  & ptr, math::Point const & p) {
                    ptr.move(p); } )
                .def("to_str",[] (motif::Motif  & ptr) -> String const {
                    return ptr.to_str(); } )
                .def("get_state",[] (motif::Motif  & ptr) -> MotifStateOP {
                    return ptr.get_state(); } )
                .def("new_res_uuids",[] (motif::Motif  & ptr) {
                    ptr.new_res_uuids(); } )
                .def("copy_uuids_from_motif",[] (motif::Motif  & ptr, motif::Motif const & m) {
                    ptr.copy_uuids_from_motif(m); } )
                .def("sequence",[] (motif::Motif  & ptr) -> String {
                    return ptr.sequence(); } )
                .def("dot_bracket",[] (motif::Motif  & ptr) -> String {
                    return ptr.dot_bracket(); } )
                .def("mtype",[] (motif::Motif  & ptr) ->  util::MotifType const & {
                    return ptr.mtype(); } )
                .def("secondary_structure",[] (motif::Motif  & ptr) ->  secondary_structure::MotifOP const & {
                    return ptr.secondary_structure(); } )
                .def("block_end_add",[] (motif::Motif  & ptr) ->  int const & {
                    return ptr.block_end_add(); } )
                .def("id",[] (motif::Motif  & ptr) ->  util::Uuid const & {
                    return ptr.id(); } )
                .def("end_name",[] (motif::Motif  & ptr, int i) ->  String {
                    return ptr.end_name(i); } )
                .def("id",[] (motif::Motif  & ptr, util::Uuid const & nid) {
                    ptr.id(nid); } )
                .def("mtype",[] (motif::Motif  & ptr, util::MotifType const & mtype) {
                    ptr.mtype(mtype); } )
                .def("secondary_structure",[] (motif::Motif  & ptr, motif::MotifOP const & ss) {
                    ptr.secondary_structure(); } )
                .def("structure",[] (motif::Motif  & ptr, structure::StructureOP const & s) {
                    ptr.structure(s); } )
                .def("block_end_add",[] (motif::Motif  & ptr, int nblock_end_add) {
                    ptr.block_end_add(nblock_end_add); } )
                .def("remove_bad_bps",[] (motif::Motif  & ptr, const structure::BasepairOPs & bad_bps) {
                    ptr.remove_bad_bps(bad_bps); } )
                .def("get_basepair",[] (motif::Motif  & ptr, String const & str) -> structure::BasepairOPs {
                    return ptr.get_basepair(str); } )
                .def("get_basepair",[] (motif::Motif  & ptr, util::Uuid const & id) -> structure::BasepairOPs {
                    return ptr.get_basepair(id); } )
                .def("get_basepair",[] (motif::Motif  & ptr, structure::ResidueOP const & r1, structure::ResidueOP const & r2) -> structure::BasepairOPs {
                    return ptr.get_basepair(r1, r2); } )
                .def("get_basepair",[] (motif::Motif  & ptr, util::Uuid const & id1, util::Uuid const & id2) -> structure::BasepairOPs {
                    return ptr.get_basepair(id1, id2); } )
                .def("get_beads",[] (motif::Motif  & ptr, structure::BasepairOPs const & bps) -> structure::Beads const & {
                    return ptr.get_beads(bps); } )
                .def("get_beads",[] (motif::Motif  & ptr, structure::BasepairOP const & bps) -> structure::Beads const & {
                    return ptr.get_beads(bps); } )
                .def("get_beads",[] (motif::Motif  & ptr) ->  structure::Beads const & {
                    return ptr.get_beads(); } )
                .def("get_end_index",[] (motif::Motif  & ptr, structure::BasepairOP const & bp) -> int {
                    return ptr.get_end_index(bp); } )
                .def("get_end_index",[] (motif::Motif  & ptr, String const & str) -> int {
                    return ptr.get_end_index(str); } )
                .def("to_pdb_str",[] (motif::Motif  & ptr, int rnumber, int close_chains) -> String const {
                    return ptr.to_pdb_str(rnumber, close_chains); } )
                .def("to_pdb",[] (motif::Motif  & ptr, String const fname, int renumber, int close_chains, int conect_statements) {
                    ptr.to_pdb(fname, renumber, close_chains, conect_statements); } )
                .def("atoms",[] (motif::Motif const & ptr) ->  structure::AtomOPs const {
                    return ptr.atoms(); } )
                .def("residues",[] (motif::Motif const & ptr) ->  structure::ResidueOPs const {
                    return ptr.residues(); } )
                .def("chains",[] (motif::Motif  & ptr) -> structure::ChainOPs const & {
                    return ptr.chains(); } )
                .def("get_residue",[] (motif::Motif  & ptr, int num, String const & chain_id, String const & i_code) ->  structure::ResidueOP const {
                    return ptr.get_residue(num, chain_id, i_code); } )
                .def("get_residue",[] (motif::Motif  & ptr, util::Uuid const & uuid) ->  structure::ResidueOP const {
                    return ptr.get_residue(uuid); } )
                .def("ends",[] (motif::Motif const & ptr) ->  structure::BasepairOPs const & {
                    return ptr.ends(); } )
                .def("end_ids",[] (motif::Motif const & ptr) ->  Strings const & {
                    return ptr.end_ids(); } )
                .def("name",[] (motif::Motif const & ptr) ->  String const & {
                    return ptr.name(); } )
                .def("path",[] (motif::Motif const & ptr) ->  String const & {
                    return ptr.path(); } )
                .def("basepairs",[] (motif::Motif const & ptr) ->  structure::BasepairOPs const & {
                    return ptr.basepairs(); } )
                .def("UNSAFE_basepairs",[] (motif::Motif  & ptr) ->  structure::BasepairOPs & {
                    return ptr.UNSAFE_basepairs(); } )
                .def("beads",[] (motif::Motif const & ptr) ->  structure::Beads const & {
                    return ptr.beads(); } )
                .def("score",[] (motif::Motif const & ptr) ->  float const & {
                    return ptr.score(); } )
                .def("protein_beads",[] (motif::Motif  & ptr) ->  structure::Beads const & {
                    return ptr.protein_beads(); } )
                .def("name",[] (motif::Motif  & ptr, String const & nname) {
                    ptr.name(nname); } )
                .def("path",[] (motif::Motif  & ptr, String const & npath) {
                    ptr.path(npath); } )
                .def("score",[] (motif::Motif  & ptr, float const & nscore) {
                    ptr.score(nscore); } )
                .def("end_ids",[] (motif::Motif  & ptr, Strings const & end_ids) {
                    ptr.end_ids(end_ids); } )
                .def("ends",[] (motif::Motif  & ptr, structure::BasepairOPs const & ends) {
                    ptr.ends(ends); } )
                .def("protein_beads",[] (motif::Motif  & ptr, structure::Beads const & beads) {
                    ptr.protein_beads(beads); } )

                ;

        py::class_<MotifEnsemble, std::shared_ptr<MotifEnsemble>>(m, "MotifEnsemble")
                // ctors
                .def(py::init<>())
                .def(py::init<motif::MotifOPs const &,Floats const &>())
                .def(py::init<String const &,motif::MotifOPs const &,Floats const &>())
                .def(py::init<MotifEnsemble const &>())
                .def(py::init<String const &,structure::ResidueTypeSet const &>())
                        // methods
                .def("size",[] (MotifEnsemble  & ptr) -> size_t {
                    return ptr.size(); } )
                .def("to_str",[] (MotifEnsemble  & ptr) -> String {
                    return ptr.to_str(); } )
                .def("get_state",[] (MotifEnsemble  & ptr) -> MotifStateEnsembleOP {
                    return ptr.get_state(); } )
                .def("members",[] (MotifEnsemble  & ptr)  {
                    return ptr.members(); } )
                .def("id",[] (MotifEnsemble  & ptr, String const & id) {
                    ptr.id(id); } )
                .def("id",[] (MotifEnsemble const & ptr) -> String {
                    return ptr.id(); } )
                ;
/*
        py::class_<MotifEnsembleMember, std::shared_ptr<MotifEnsembleMember>>(m, "MotifEnsembleMember")
		// ctors
		.def(py::init<MotifOP const &,float const &>())
		// methods
		.def("to_str",[] (MotifEnsembleMember  & ptr) -> String {
		 return ptr.to_str(); } )
		// public attributes
		.def_readwrite("motif", &MotifEnsembleMember::motif)
		.def_readwrite("energy", &MotifEnsembleMember::energy)
		;
*/
        py::class_<MotifFactory, std::shared_ptr<MotifFactory>>(m, "MotifFactory")
                // ctors
                .def(py::init<>())
                        // methods
                .def("motif_from_file",[] (MotifFactory  & ptr, String const & path, bool rebuild_x3dna, bool include_protein, int force_num_chains)  {
                    return *ptr.motif_from_file(path, rebuild_x3dna, include_protein, force_num_chains); } )
                .def("motif_from_res",[] (MotifFactory  & ptr, structure::ResidueOPs & res, structure::BasepairOPs const & bps)  {
                    return *ptr.motif_from_res(res, bps); } )
                .def("motif_from_bps",[] (MotifFactory  & ptr, structure::BasepairOPs const & bps){
                    return *ptr.motif_from_bps(bps); } )
                .def("can_align_motif_to_end",[] (MotifFactory  & ptr, motif::MotifOP const & m, int ei)  {
                    return *ptr.can_align_motif_to_end(m, ei); } )
                .def("align_motif_to_common_frame",[] (MotifFactory  & ptr, motif::MotifOP const & m, int ei){
                    return *ptr.align_motif_to_common_frame(m, ei); } )
                .def("standardize_rna_structure_ends",[] (MotifFactory  & ptr, motif::MotifOP & m) {
                    ptr.standardize_rna_structure_ends(m); } )
                .def("_setup_basepair_ends",[] (MotifFactory  & ptr, structure::StructureOP const & structure, structure::BasepairOPs const & basepairs) -> structure::BasepairOPs {
                    return ptr._setup_basepair_ends(structure, basepairs); } )
                .def("_setup_secondary_structure",[] (MotifFactory  & ptr, motif::MotifOP & m) {
                    ptr._setup_secondary_structure(m); } )
//            .def("added_helix",[] (MotifFactory  & ptr)  {
//                return *ptr.added_helix(); }
//                )
                .def("ref_motif",[] (MotifFactory  & ptr) {
                    return *ptr.ref_motif(); } )
                ;

        py::class_<MotifScorer, std::shared_ptr<MotifScorer>>(m, "MotifScorer")
                // ctors
                .def(py::init<>())
                        // methods
                .def("score",[] (MotifScorer  & ptr, motif::MotifOP const & m) -> float  {
                    return ptr.score(m); } )
                ;

        py::class_<MotifState, std::shared_ptr<MotifState>>(m, "MotifState")
                // ctors
                .def(py::init<String const &,Strings const &,Strings const &,structure::BasepairStateOPs const &,math::Points const &,float const &,int const &,int const &>())
                .def(py::init<String const &,Strings const &,Strings const &,structure::BasepairStateOPs const &,math::Points const &,float const &,int const &,int const &,util::Uuid const &>())
                .def(py::init<MotifState const &>())
                .def(py::init<String const &>())
                        // methods
                .def("move",[] (MotifState  & ptr, math::Point const & p) {
                    ptr.move(p); } )
                .def("transform",[] (MotifState  & ptr, math::Matrix const & r, math::Vector const & t, math::Point & dummy) {
                    ptr.transform(r, t, dummy); } )
                .def("transform",[] (MotifState  & ptr, math::Matrix const & r, math::Vector const & t) {
                    ptr.transform(r, t); } )
                .def("to_str",[] (MotifState  & ptr) -> String {
                    return ptr.to_str(); } )
                .def("get_end_state",[] (MotifState  & ptr, String const & name)  {
                    return ptr.get_end_state(name); } )
                .def("get_end_index",[] (MotifState  & ptr, structure::BasepairStateOP const & end_state) -> int {
                    return ptr.get_end_index(end_state); } )
                .def("get_end_index",[] (MotifState  & ptr, structure::BasepairStateOP const & end_state) -> int  {
                    return ptr.get_end_index(end_state); } )
                .def("name",[] (MotifState  & ptr) ->  String const & {
                    return ptr.name(); } )
                .def("end_names",[] (MotifState const & ptr) ->  Strings const & {
                    return ptr.end_names(); } )
                .def("end_ids",[] (MotifState  & ptr) ->  Strings const & {
                    return ptr.end_ids(); } )
                .def("end_states",[] (MotifState  & ptr) ->  structure::BasepairStateOPs & {
                    return ptr.end_states(); } )
                .def("beads",[] (MotifState const & ptr) ->  math::Points const & {
                    return ptr.beads(); } )
                .def("score",[] (MotifState  & ptr) ->  float const & {
                    return ptr.score(); } )
                .def("size",[] (MotifState  & ptr) ->  int const & {
                    return ptr.size(); } )
                .def("block_end_add",[] (MotifState  & ptr) ->  int const & {
                    return ptr.block_end_add(); } )
                .def("uuid",[] (MotifState  & ptr) ->  util::Uuid const & {
                    return ptr.uuid(); } )
                .def("uuid",[] (MotifState  & ptr, util::Uuid const & uuid) {
                    ptr.uuid(uuid); } )
                .def("update_end_state",[] (MotifState  & ptr, int i, structure::BasepairState const & new_state) {
                    ptr.update_end_state(i, new_state); } )
                .def("beads",[] (MotifState  & ptr, math::Points const & beads) {
                    ptr.beads(beads); } )
                .def("new_uuids",[] (MotifState  & ptr) {
                    ptr.new_uuids(); } )
                ;

        py::class_<MotifStateAligner, std::shared_ptr<MotifStateAligner>>(m, "MotifStateAligner")
                // ctors
                .def(py::init<>())
                        // methods
                .def("get_aligned_motif_state",[] (MotifStateAligner  & ptr, structure::BasepairStateOP const & ref_bp_state, MotifStateOP & cur_state, MotifStateOP const & org_state) {
                    ptr.get_aligned_motif_state(ref_bp_state, cur_state, org_state); } )
                .def("get_aligned_motif_state",[] (MotifStateAligner  & ptr, structure::BasepairStateOP ref_bp_state, MotifStateOP state) {
                    ptr.get_aligned_motif_state(ref_bp_state, state); } )
                ;

        py::class_<MotifStateEnsemble, std::shared_ptr<MotifStateEnsemble>>(m, "MotifStateEnsemble")
                // ctors
                .def(py::init<MotifStateOP const &>())
                .def(py::init<MotifStateEnsemble const &>())
                .def(py::init<MotifStateOPs const,Floats const>())
                .def(py::init<String const>())
                        // methods
                .def("begin",[] (MotifStateEnsemble  & ptr) -> MotifStateEnsembleMemberOPs::iterator {
                    return ptr.begin(); } )
                .def("end",[] (MotifStateEnsemble  & ptr) -> MotifStateEnsembleMemberOPs::iterator {
                    return ptr.end(); } )
                .def("begin",[] (MotifStateEnsemble const & ptr) -> MotifStateEnsembleMemberOPs::const_iterator {
                    return ptr.begin(); } )
                .def("end",[] (MotifStateEnsemble const & ptr) -> MotifStateEnsembleMemberOPs::const_iterator {
                    return ptr.end(); } )
                .def("to_str",[] (MotifStateEnsemble  & ptr) -> String {
                    return ptr.to_str(); } )
                .def("num_end_states",[] (MotifStateEnsemble const & ptr) ->  int {
                    return ptr.num_end_states(); } )
                .def("most_populated",[] (MotifStateEnsemble const & ptr) ->  MotifStateOP const & {
                    return ptr.most_populated(); } )
                .def("get_member",[] (MotifStateEnsemble const & ptr, int i) ->  MotifStateEnsembleMemberOP const & {
                    return ptr.get_member(i); } )
                .def("member_index",[] (MotifStateEnsemble  & ptr, motif::MotifStateEnsembleMemberOP mem) ->  int {
                    return ptr.member_index(mem); } )
                .def("members",[] (MotifStateEnsemble  & ptr) ->  MotifStateEnsembleMemberOPs const & {
                    return ptr.members(); } )
                .def("size",[] (MotifStateEnsemble const & ptr) -> size_t {
                    return ptr.size(); } )
                .def("block_end_add",[] (MotifStateEnsemble  & ptr) ->  int {
                    return ptr.block_end_add(); } )
                .def("id",[] (MotifStateEnsemble  & ptr) ->  String {
                    return ptr.id(); } )
                ;

        py::class_<MotifStateEnsembleMember, std::shared_ptr<MotifStateEnsembleMember>>(m, "MotifStateEnsembleMember")
                // ctors
                .def(py::init<MotifStateOP const &,float const &>())
                .def(py::init<MotifStateEnsembleMember const &>())
                        // methods
                .def("to_str",[] (MotifStateEnsembleMember  & ptr) ->  String {
                    return ptr.to_str(); } )
                        // public attributes
                .def_readwrite("motif_state", &MotifStateEnsembleMember::motif_state)
                .def_readwrite("energy", &MotifStateEnsembleMember::energy)
                ;
/*
        py::class_<MotifStateEnsembleMember_LessThanKey, std::shared_ptr<MotifStateEnsembleMember_LessThanKey>>(m, "MotifStateEnsembleMember_LessThanKey")
		// operators
		.def(py::self () motif::MotifStateEnsembleMemberOP)
		;
*/
        py::class_<MotiftoSecondaryStructure, std::shared_ptr<MotiftoSecondaryStructure>>(m, "MotiftoSecondaryStructure")
                // ctors
                .def(py::init<>())
                        // methods
                .def("reset",[] (MotiftoSecondaryStructure  & ptr) {
                    ptr.reset(); } )
                .def("to_secondary_structure",[] (MotiftoSecondaryStructure  & ptr, structure::RNAStructureOP const & motif) -> secondary_structure::RNAStructureOP {
                    return ptr.to_secondary_structure(motif); } )
                ;

        py::class_<motif::Pose, std::shared_ptr<motif::Pose>>(m, "Pose")
                // ctors
                .def(py::init<>())
                .def(py::init<motif::MotifOP const>())
                .def(py::init<structure::StructureOP const,structure::BasepairOPs const>())
                        // methods
                .def("motifs",[] (motif::Pose  & ptr, util::MotifType const & mtype) -> motif::MotifOPs const & {
                    return ptr.motifs(mtype); } )
                .def("set_bp_designable",[] (motif::Pose  & ptr, structure::BasepairOP const & bp) {
                    ptr.set_bp_designable(bp); } )
                .def("designable_sequence",[] (motif::Pose  & ptr) -> String {
                    return ptr.designable_sequence(); } )
                .def("designable",[] (motif::Pose  & ptr, std::map<util::Uuid, int, util::UuidCompare> const & ndesignable) {
                    ptr.designable(ndesignable); } )
                .def("set_motifs",[] (motif::Pose  & ptr, std::map<util::MotifType, motif::MotifOPs> const & motifs) {
                    ptr.set_motifs(motifs); } )
                        // inherited methods
                .def("transform",[] (motif::Motif  & ptr, math::Transform const & t) {
                    ptr.transform(t); } )
                .def("move",[] (motif::Motif  & ptr, math::Point const & p) {
                    ptr.move(p); } )
                .def("to_str",[] (motif::Motif  & ptr) -> String const {
                    return ptr.to_str(); } )
                .def("get_state",[] (motif::Motif  & ptr) -> motif::MotifStateOP {
                    return ptr.get_state(); } )
                .def("new_res_uuids",[] (motif::Motif  & ptr) {
                    ptr.new_res_uuids(); } )
                .def("copy_uuids_from_motif",[] (motif::Motif  & ptr, motif::Motif const & m) {
                    ptr.copy_uuids_from_motif(m); } )
                .def("sequence",[] (motif::Motif  & ptr) -> String {
                    return ptr.sequence(); } )
                .def("dot_bracket",[] (motif::Motif  & ptr) -> String {
                    return ptr.dot_bracket(); } )
                .def("mtype",[] (motif::Motif  & ptr) ->  util::MotifType const & {
                    return ptr.mtype(); } )
                .def("secondary_structure",[] (motif::Motif  & ptr) ->  secondary_structure::MotifOP const & {
                    return ptr.secondary_structure(); } )
                .def("block_end_add",[] (motif::Motif  & ptr) ->  int const & {
                    return ptr.block_end_add(); } )
                .def("id",[] (motif::Motif  & ptr) ->  util::Uuid const & {
                    return ptr.id(); } )
                .def("end_name",[] (motif::Motif  & ptr, int i) ->  String {
                    return ptr.end_name(i); } )
                .def("id",[] (motif::Motif  & ptr, util::Uuid const & nid) {
                    ptr.id(nid); } )
                .def("mtype",[] (motif::Motif  & ptr, util::MotifType const & mtype) {
                    ptr.mtype(mtype); } )
                .def("secondary_structure",[] (motif::Motif  & ptr, secondary_structure::MotifOP const & ss) {
                    ptr.secondary_structure(ss); } )
                .def("structure",[] (motif::Motif  & ptr, structure::StructureOP const & s) {
                    ptr.structure(s); } )
                .def("block_end_add",[] (motif::Motif  & ptr, int nblock_end_add) {
                    ptr.block_end_add(nblock_end_add); } )
                .def("remove_bad_bps",[] (motif::Motif  & ptr, const structure::BasepairOPs & bad_bps) {
                    ptr.remove_bad_bps(bad_bps); } )
                .def("get_basepair",[] (motif::Pose  & ptr, String const & str) {
                    return ptr.get_basepair(str); } )
                .def("get_basepair",[] (motif::Pose  & ptr, util::Uuid const & id){
                    return ptr.get_basepair(id); } )
                .def("get_basepair",[] (motif::Pose  & ptr, structure::ResidueOP const & r1, structure::ResidueOP const & r2) {
                    return ptr.get_basepair(r1, r2); } )
                .def("get_basepair",[] (motif::Pose  & ptr, util::Uuid const & id1, util::Uuid const & id2)  {
                    return ptr.get_basepair(id1, id2); } )
                .def("get_beads",[] (motif::Pose  & ptr, structure::BasepairOPs const & bps) -> structure::Beads const & {
                    return ptr.get_beads(bps); } )
                .def("get_beads",[] (motif::Pose  & ptr, structure::BasepairOP const & bps) -> structure::Beads const & {
                    return ptr.get_beads(bps); } )
                .def("get_beads",[] (motif::Pose  & ptr) ->  structure::Beads const & {
                    return ptr.get_beads(); } )
                .def("get_end_index",[] (motif::Pose  & ptr, structure::BasepairOP const & bp) -> int {
                    return ptr.get_end_index(bp); } )
                .def("get_end_index",[] (motif::Pose  & ptr, String const & str) -> int {
                    return ptr.get_end_index(str); } )
                .def("to_pdb_str",[] (motif::Pose  & ptr, int rnumber, int close_chains) -> String const {
                    return ptr.to_pdb_str(rnumber, close_chains); } )
                .def("to_pdb",[] (motif::Pose  & ptr, String const fname, int renumber, int close_chains, int conect_statements) {
                    ptr.to_pdb(fname, renumber, close_chains, conect_statements); } )
                .def("atoms",[] (motif::Pose const & ptr) ->  structure::AtomOPs const {
                    return ptr.atoms(); } )
                .def("residues",[] (motif::Pose const & ptr) ->  structure::ResidueOPs const {
                    return ptr.residues(); } )
                .def("chains",[] (motif::Pose  & ptr) -> structure::ChainOPs const & {
                    return ptr.chains(); } )
                .def("get_residue",[] (motif::Pose  & ptr, int num, String const & chain_id, String const & i_code) ->  structure::ResidueOP const {
                    return ptr.get_residue(num, chain_id, i_code); } )
                .def("get_residue",[] (motif::Pose  & ptr, util::Uuid const & uuid) ->  structure::ResidueOP const {
                    return ptr.get_residue(uuid); } )
                .def("ends",[] (motif::Pose const & ptr) ->  structure::BasepairOPs const & {
                    return ptr.ends(); } )
                .def("end_ids",[] (motif::Pose const & ptr) ->  Strings const & {
                    return ptr.end_ids(); } )
                .def("name",[] (motif::Pose const & ptr) ->  String const & {
                    return ptr.name(); } )
                .def("path",[] (motif::Pose const & ptr) ->  String const & {
                    return ptr.path(); } )
                .def("basepairs",[] (motif::Pose const & ptr) ->  structure::BasepairOPs const & {
                    return ptr.basepairs(); } )
                .def("UNSAFE_basepairs",[] (motif::Pose  & ptr) ->  structure::BasepairOPs & {
                    return ptr.UNSAFE_basepairs(); } )
                .def("beads",[] (motif::Pose const & ptr) ->  structure::Beads const & {
                    return ptr.beads(); } )
                .def("score",[] (motif::Pose const & ptr) ->  float const & {
                    return ptr.score(); } )
                .def("protein_beads",[] (motif::Pose  & ptr) ->  structure::Beads const & {
                    return ptr.protein_beads(); } )
                .def("name",[] (motif::Pose  & ptr, String const & nname) {
                    ptr.name(nname); } )
                .def("path",[] (motif::Pose  & ptr, String const & npath) {
                    ptr.path(npath); } )
                .def("score",[] (motif::Pose  & ptr, float const & nscore) {
                    ptr.score(nscore); } )
                .def("end_ids",[] (motif::Pose  & ptr, Strings const & end_ids) {
                    ptr.end_ids(end_ids); } )
                .def("ends",[] (motif::Pose  & ptr, structure::BasepairOPs const & ends) {
                    ptr.ends(ends); } )
                .def("protein_beads",[] (motif::Pose  & ptr, structure::Beads const & beads) {
                    ptr.protein_beads(beads); } )

                ;

        py::class_<PoseFactory, std::shared_ptr<PoseFactory>>(m, "PoseFactory")
                // ctors
                .def(py::init<>())
                        // methods
                .def("pose_from_motif_tree",[] (PoseFactory  & ptr, structure::StructureOP const & structure, structure::BasepairOPs const & basepairs, motif::MotifOPs const & motifs, std::map<util::Uuid, int, util::UuidCompare> const & designable) -> motif::PoseOP {
                    return ptr.pose_from_motif_tree(structure, basepairs, motifs, designable); } )
                .def("pose_from_file",[] (PoseFactory  & ptr, String const & path, int gu_are_helix, int signlet_bp_seperation) -> motif::PoseOP {
                    return ptr.pose_from_file(path, gu_are_helix, signlet_bp_seperation); } )
                ;

    }

}

#endif // PYBIND11_MOTIF_HPP