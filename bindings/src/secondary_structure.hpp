#ifndef PYBIND11_SECONDARY_STRUCTURE_HPP
#define PYBIND11_SECONDARY_STRUCTURE_HPP

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/cast.h>
#include <pybind11/operators.h>
#include <memory>

// secondary_structure
#include <secondary_structure/basepair.h>
#include <secondary_structure/chain.h>
#include <secondary_structure/motif.h>
#include <secondary_structure/pose.h>
#include <secondary_structure/residue.h>
#include <secondary_structure/rna_structure.h>
#include <secondary_structure/secondary_structure_parser.h>
#include <secondary_structure/secondary_structure_tree.h>
#include <secondary_structure/sequence_constraint.h>
#include <secondary_structure/sequence_tools.h>
#include <secondary_structure/structure.h>
#include <secondary_structure/util.h>

namespace secondary_structure {
    namespace py = pybind11;
    void
    add_bindings(py::module_& m) {
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// secondary_structure
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // free functions
        m.def("assign_end_id", [] (secondary_structure::RNAStructureOP const & ss, secondary_structure::BasepairOP const & end) -> String {
                          return assign_end_id(ss, end); },
                      py::arg("ss"),
                      py::arg("end")
        );

        m.def("convert_res_name_to_type", [] (char c) -> ResType {
                          return convert_res_name_to_type(c); },
                      py::arg("c")
        );

        m.def("fill_basepairs_in_ss", [] (secondary_structure::PoseOP & ss) {
                          fill_basepairs_in_ss(ss); },
                      py::arg("ss")
        );

        m.def("find_gc_helix_stretches", [] (secondary_structure::PoseOP p, int length) -> int {
                          return find_gc_helix_stretches(p, length); },
                      py::arg("p"),
                      py::arg("length")
        );

        m.def("find_longest_gc_helix_stretch", [] (secondary_structure::PoseOP p) -> int {
                          return find_longest_gc_helix_stretch(p); },
                      py::arg("p")
        );

        m.def("find_res_types_in_pose", [] (secondary_structure::PoseOP p, ResTypes const & residue_types) -> int {
                          return find_res_types_in_pose(p, residue_types); },
                      py::arg("p"),
                      py::arg("residue_types")
        );

        m.def("get_res_types_from_sequence", [] (String const & sequence, ResTypes & residue_types) {
                          get_res_types_from_sequence(sequence, residue_types); },
                      py::arg("sequence"),
                      py::arg("residue_types")
        );

        m.def("is_au_pair", [] (secondary_structure::BasepairOP bp) ->  bool {
                  return is_au_pair(bp); },
              py::arg("bp")
        );

        m.def("is_gc_pair", [] (secondary_structure::BasepairOP bp) ->  bool {
                  return is_gc_pair(bp); },
              py::arg("bp")
        );

        m.def("is_gu_pair", [] (secondary_structure::BasepairOP bp) ->  bool {
                  return is_gu_pair(bp); },
              py::arg("bp")
        );

        m.def("tree_from_pose", [] (secondary_structure::PoseOP const & p) -> SecondaryStructureTreeOP {
                  return tree_from_pose(p); },
              py::arg("p")
        );
        // classes

        py::class_<secondary_structure::Basepair, std::shared_ptr<secondary_structure::Basepair>>(m, "Basepair")
                // ctors
                .def(py::init<secondary_structure::ResidueOP const &,secondary_structure::ResidueOP const &,util::Uuid const &>())
                        // methods
                .def("name",[] (secondary_structure::Basepair  & ptr) ->  String {
                    return ptr.name(); } )
                .def("partner",[] (secondary_structure::Basepair  & ptr, secondary_structure::ResidueOP const & r) ->  secondary_structure::ResidueOP {
                    return ptr.partner(r); } )
                .def("res1",[] (secondary_structure::Basepair  & ptr) ->  secondary_structure::ResidueOP & {
                    return ptr.res1(); } )
                .def("res2",[] (secondary_structure::Basepair  & ptr) ->  secondary_structure::ResidueOP & {
                    return ptr.res2(); } )
                .def("res1",[] (secondary_structure::Basepair const & ptr) ->  secondary_structure::ResidueOP const & {
                    return ptr.res1(); } )
                .def("res2",[] (secondary_structure::Basepair const & ptr) ->  secondary_structure::ResidueOP const & {
                    return ptr.res2(); } )
                .def("uuid",[] (secondary_structure::Basepair  & ptr) ->  util::Uuid const & {
                    return ptr.uuid(); } )
                ;

        py::class_<secondary_structure::Chain, std::shared_ptr<secondary_structure::Chain>>(m, "Chain")
                // ctors
                .def(py::init<>())
                .def(py::init<secondary_structure::ResidueOPs const &>())
                .def(py::init<secondary_structure::Chain const &>())
                .def(py::init<String const &>())
                        // methods
                .def("begin",[] (secondary_structure::Chain  & ptr) -> secondary_structure::ResidueOPs::iterator {
                    return ptr.begin(); } )
                .def("end",[] (secondary_structure::Chain  & ptr) -> secondary_structure::ResidueOPs::iterator {
                    return ptr.end(); } )
                .def("begin",[] (secondary_structure::Chain const & ptr) -> secondary_structure::ResidueOPs::const_iterator {
                    return ptr.begin(); } )
                .def("end",[] (secondary_structure::Chain const & ptr) -> secondary_structure::ResidueOPs::const_iterator {
                    return ptr.end(); } )
                .def("first",[] (secondary_structure::Chain  & ptr) ->  secondary_structure::ResidueOP const & {
                    return ptr.first(); } )
                .def("last",[] (secondary_structure::Chain  & ptr) ->  secondary_structure::ResidueOP const & {
                    return ptr.last(); } )
                .def("sequence",[] (secondary_structure::Chain  & ptr) ->  String {
                    return ptr.sequence(); } )
                .def("dot_bracket",[] (secondary_structure::Chain  & ptr) ->  String {
                    return ptr.dot_bracket(); } )
                .def("to_str",[] (secondary_structure::Chain  & ptr) ->  String {
                    return ptr.to_str(); } )
                .def("length",[] (secondary_structure::Chain  & ptr) ->  int {
                    return ptr.length(); } )
                .def("residues",[] (secondary_structure::Chain  & ptr) ->  secondary_structure::ResidueOPs const & {
                    return ptr.residues(); } )
                ;

        py::class_<DisallowedSequence, std::shared_ptr<DisallowedSequence>>(m, "DisallowedSequence")
                // ctors
                .def(py::init<String const &>())
                        // methods
                .def("clone",[] (DisallowedSequence const & ptr) -> SequenceConstraint * {
                    return ptr.clone(); } )
                .def("violations",[] (DisallowedSequence  & ptr, secondary_structure::PoseOP p) -> int {
                    return ptr.violations(p); } )
                        // inherited methods
                .def("clone",[] (SequenceConstraint const & ptr) -> SequenceConstraint * {
                    return ptr.clone(); } )
                .def("violates_constraint",[] (SequenceConstraint  & ptr, secondary_structure::PoseOP p) -> bool {
                    return ptr.violates_constraint(p); } )
                .def("violations",[] (SequenceConstraint  & ptr, secondary_structure::PoseOP p) -> int {
                    return ptr.violations(p); } )
                ;

        py::class_<GCHelixStretchLimit, std::shared_ptr<GCHelixStretchLimit>>(m, "GCHelixStretchLimit")
                // ctors
                .def(py::init<int>())
                        // methods
                .def("clone",[] (GCHelixStretchLimit const & ptr) -> SequenceConstraint * {
                    return ptr.clone(); } )
                .def("violations",[] (GCHelixStretchLimit  & ptr, secondary_structure::PoseOP p) -> int {
                    return ptr.violations(p); } )
                        // inherited methods
                .def("clone",[] (SequenceConstraint const & ptr) -> SequenceConstraint * {
                    return ptr.clone(); } )
                .def("violates_constraint",[] (SequenceConstraint  & ptr, secondary_structure::PoseOP p) -> bool {
                    return ptr.violates_constraint(p); } )
                .def("violations",[] (SequenceConstraint  & ptr, secondary_structure::PoseOP p) -> int {
                    return ptr.violations(p); } )
                ;

        py::class_<secondary_structure::Motif, std::shared_ptr<secondary_structure::Motif>>(m, "Motif")
                // ctors
                .def(py::init<>())
                .def(py::init<secondary_structure::StructureOP const &,secondary_structure:: BasepairOPs const &,secondary_structure::BasepairOPs const &>())
                .def(py::init<secondary_structure::StructureOP const &,secondary_structure:: BasepairOPs const &,secondary_structure::BasepairOPs const &,Strings const &,String const &,String const &,float>())
                .def(py::init<secondary_structure::Motif const &>())
                .def(py::init<String const &>())
                        // methods
                .def("to_str",[] (secondary_structure::Motif  & ptr) -> String {
                    return ptr.to_str(); } )
                .def("mtype",[] (secondary_structure::Motif  & ptr) ->  util::MotifType const & {
                    return ptr.mtype(); } )
                .def("id",[] (secondary_structure::Motif  & ptr) ->  util::Uuid const & {
                    return ptr.id(); } )
                .def("mtype",[] (secondary_structure::Motif  & ptr, util::MotifType const & mtype) {
                    ptr.mtype(mtype); } )
                .def("id",[] (secondary_structure::Motif  & ptr, util::Uuid const & uuid) {
                    ptr.id(uuid); } )
                        // inherited methods
                .def("get_basepair",[] (secondary_structure::RNAStructure  & ptr, util::Uuid const & bp_uuid) ->secondary_structure::BasepairOPs {
                    return ptr.get_basepair(bp_uuid); } )
                .def("get_basepair",[] (secondary_structure::RNAStructure  & ptr, util::Uuid const & bp_uuid) ->secondary_structure::BasepairOPs {
                    return ptr.get_basepair(bp_uuid); } )
                .def("get_basepair",[] (secondary_structure::RNAStructure  & ptr, util::Uuid const & bp_uuid) ->secondary_structure::BasepairOPs {
                    return ptr.get_basepair(bp_uuid); } )
                .def("get_basepair",[] (secondary_structure::RNAStructure  & ptr, util::Uuid const & bp_uuid) ->secondary_structure::BasepairOPs {
                    return ptr.get_basepair(bp_uuid); } )
                .def("get_end",[] (secondary_structure::RNAStructure  & ptr, String const & name) -> secondary_structure::BasepairOP {
                    return ptr.get_end(name); } )
                .def("replace_sequence",[] (secondary_structure::RNAStructure  & ptr, String const & seq) {
                    ptr.replace_sequence(seq); } )
                .def("get_residue",[] (secondary_structure::RNAStructure  & ptr, int num, String const & chain_id, String const & i_code) ->  secondary_structure::ResidueOP {
                    return ptr.get_residue(num, chain_id, i_code); } )
                .def("get_residue",[] (secondary_structure::RNAStructure  & ptr, util::Uuid const & uuid) ->  secondary_structure::ResidueOP {
                    return ptr.get_residue(uuid); } )
                .def("sequence",[] (secondary_structure::RNAStructure  & ptr) ->  String {
                    return ptr.sequence(); } )
                .def("dot_bracket",[] (secondary_structure::RNAStructure  & ptr) ->  String {
                    return ptr.dot_bracket(); } )
                .def("chains",[] (secondary_structure::RNAStructure  & ptr) ->  secondary_structure::ChainOPs const & {
                    return ptr.chains(); } )
                .def("residues",[] (secondary_structure::RNAStructure  & ptr) ->  secondary_structure::ResidueOPs {
                    return ptr.residues(); } )
                .def("structure",[] (secondary_structure::RNAStructure  & ptr) ->  secondary_structure::StructureOP {
                    return ptr.structure(); } )
                .def("basepairs",[] (secondary_structure::RNAStructure  & ptr) ->  secondary_structure::BasepairOPs const & {
                    return ptr.basepairs(); } )
                .def("ends",[] (secondary_structure::RNAStructure  & ptr) ->  secondary_structure::BasepairOPs const & {
                    return ptr.ends(); } )
                .def("name",[] (secondary_structure::RNAStructure  & ptr) ->  String const & {
                    return ptr.name(); } )
                .def("end_ids",[] (secondary_structure::RNAStructure  & ptr) ->  Strings const & {
                    return ptr.end_ids(); } )
                .def("name",[] (secondary_structure::RNAStructure  & ptr, String const & name) {
                    ptr.name(name); } )
                .def("path",[] (secondary_structure::RNAStructure  & ptr, String const & path) {
                    ptr.path(path); } )
                .def("end_ids",[] (secondary_structure::RNAStructure  & ptr, Strings const & end_ids) {
                    ptr.end_ids(end_ids); } )
                ;

        py::class_<NodeData, std::shared_ptr<NodeData>>(m, "NodeData")
                // ctors
                .def(py::init<>())
                .def(py::init<secondary_structure::ResidueOPs const &,NodeType const &>())
                        // public attributes
                .def_readwrite("residues", &NodeData::residues)
                .def_readwrite("type", &NodeData::type)
                ;

        py::class_<Parser, std::shared_ptr<Parser>>(m, "Parser")
                // ctors
                .def(py::init<>())
                        // methods
                .def("parse",[] (Parser  & ptr, String const & sequence, String const & dot_bracket) -> SecondaryStructureChainGraphOP {
                    return ptr.parse(sequence, dot_bracket); } )
                .def("parse_to_motifs",[] (Parser  & ptr, String const & sequence, String const & dot_bracket)   {
                    return ptr.parse_to_motifs(sequence, dot_bracket); } )
                .def("parse_to_motif",[] (Parser  & ptr, String const & sequence, String const & dot_bracket) -> secondary_structure::MotifOP  {
                    return ptr.parse_to_motif(sequence, dot_bracket); } )
                .def("parse_to_pose",[] (Parser  & ptr, String const & sequence, String const & dot_bracket) ->secondary_structure:: PoseOP {
                    return ptr.parse_to_pose(sequence, dot_bracket); } )
                .def("reset",[] (Parser  & ptr) {
                    ptr.reset(); } )
                ;

        py::class_<secondary_structure::Pose, std::shared_ptr<secondary_structure::Pose>>(m, "Pose")
                // ctors
                .def(py::init<>())
                .def(py::init<secondary_structure::StructureOP const &,secondary_structure::BasepairOPs const &,secondary_structure::BasepairOPs const &>())
                .def(py::init<secondary_structure::StructureOP const &,secondary_structure::BasepairOPs const &,secondary_structure::BasepairOPs const &,secondary_structure::MotifOPs const &>())
                .def(py::init<secondary_structure::RNAStructureOP const &,secondary_structure::MotifOPs const &>())
                        // methods
                .def("helices",[] (secondary_structure::Pose  & ptr) -> secondary_structure::MotifOPs const & {
                    return ptr.helices(); } )
                .def("motifs",[] (secondary_structure::Pose const & ptr) -> secondary_structure::MotifOPs const & {
                    return ptr.motifs(); } )
                .def("motif",[] (secondary_structure::Pose  & ptr, util::Uuid const & uuid) -> secondary_structure::MotifOP {
                    return ptr.motif(uuid); } )
                .def("replace_sequence",[] (secondary_structure::Pose  & ptr, String const & seq) {
                    ptr.replace_sequence(seq); } )
                .def("update_motif",[] (secondary_structure::Pose  & ptr, util::Uuid const & uuid) {
                    ptr.update_motif(uuid); } )
                        // inherited methods
                .def("get_basepair",[] (secondary_structure::RNAStructure  & ptr, util::Uuid const & bp_uuid) -> secondary_structure::BasepairOPs {
                    return ptr.get_basepair(bp_uuid); } )
                .def("get_basepair",[] (secondary_structure::RNAStructure  & ptr, util::Uuid const & bp_uuid) -> secondary_structure::BasepairOPs {
                    return ptr.get_basepair(bp_uuid); } )
                .def("get_basepair",[] (secondary_structure::RNAStructure  & ptr, util::Uuid const & bp_uuid) -> secondary_structure::BasepairOPs  {
                    return ptr.get_basepair(bp_uuid); } )
                .def("get_basepair",[] (secondary_structure::RNAStructure  & ptr, util::Uuid const & bp_uuid) -> secondary_structure::BasepairOPs {
                    return ptr.get_basepair(bp_uuid); } )
                .def("get_end",[] (secondary_structure::RNAStructure  & ptr, String const & name) -> secondary_structure::BasepairOP {
                    return ptr.get_end(name); } )
                .def("replace_sequence",[] (secondary_structure::RNAStructure  & ptr, String const & seq) {
                    ptr.replace_sequence(seq); } )
                .def("get_residue",[] (secondary_structure::RNAStructure  & ptr, int num, String const & chain_id, String const & i_code) ->  secondary_structure::ResidueOP {
                    return ptr.get_residue(num, chain_id, i_code); } )
                .def("get_residue",[] (secondary_structure::RNAStructure  & ptr, util::Uuid const & uuid) ->  secondary_structure::ResidueOP {
                    return ptr.get_residue(uuid); } )
                .def("sequence",[] (secondary_structure::RNAStructure  & ptr) ->  String {
                    return ptr.sequence(); } )
                .def("dot_bracket",[] (secondary_structure::RNAStructure  & ptr) ->  String {
                    return ptr.dot_bracket(); } )
                .def("chains",[] (secondary_structure::RNAStructure  & ptr) ->  secondary_structure::ChainOPs const & {
                    return ptr.chains(); } )
                .def("residues",[] (secondary_structure::RNAStructure  & ptr) ->  secondary_structure::ResidueOPs {
                    return ptr.residues(); } )
                .def("structure",[] (secondary_structure::RNAStructure  & ptr) ->  secondary_structure::StructureOP {
                    return ptr.structure(); } )
                .def("basepairs",[] (secondary_structure::RNAStructure  & ptr) ->  secondary_structure::BasepairOPs const & {
                    return ptr.basepairs(); } )
                .def("ends",[] (secondary_structure::RNAStructure  & ptr) ->  secondary_structure::BasepairOPs const & {
                    return ptr.ends(); } )
                .def("name",[] (secondary_structure::RNAStructure  & ptr) ->  String const & {
                    return ptr.name(); } )
                .def("end_ids",[] (secondary_structure::RNAStructure  & ptr) ->  Strings const & {
                    return ptr.end_ids(); } )
                .def("name",[] (secondary_structure::RNAStructure  & ptr, String const & name) {
                    ptr.name(name); } )
                .def("path",[] (secondary_structure::RNAStructure  & ptr, String const & path) {
                    ptr.path(path); } )
                .def("end_ids",[] (secondary_structure::RNAStructure  & ptr, Strings const & end_ids) {
                    ptr.end_ids(end_ids); } )
                ;

        py::class_<secondary_structure::RNAStructure, std::shared_ptr<secondary_structure::RNAStructure>>(m, "RNAStructure")
                // ctors
                .def(py::init<>())
                .def(py::init<secondary_structure::StructureOP const &,secondary_structure::BasepairOPs const &,secondary_structure:: BasepairOPs const &>())
                .def(py::init<secondary_structure::StructureOP const &,secondary_structure::BasepairOPs const &,secondary_structure:: BasepairOPs const &,Strings const &,String const &,String const &,float>())
                .def(py::init<secondary_structure::RNAStructure const &>())
                        // methods
                .def("get_basepair",[] (secondary_structure::RNAStructure  & ptr, util::Uuid const & bp_uuid) -> secondary_structure::BasepairOPs {
                    return ptr.get_basepair(bp_uuid); } )
                .def("get_basepair",[] (secondary_structure::RNAStructure  & ptr, util::Uuid const & bp_uuid) -> secondary_structure::BasepairOPs {
                    return ptr.get_basepair(bp_uuid); } )
                .def("get_basepair",[] (secondary_structure::RNAStructure  & ptr, util::Uuid const & bp_uuid) -> secondary_structure::BasepairOPs {
                    return ptr.get_basepair(bp_uuid); } )
                .def("get_basepair",[] (secondary_structure::RNAStructure  & ptr, util::Uuid const & bp_uuid) -> secondary_structure::BasepairOPs {
                    return ptr.get_basepair(bp_uuid); } )
                .def("get_end",[] (secondary_structure::RNAStructure  & ptr, String const & name) -> secondary_structure::BasepairOP  {
                    return ptr.get_end(name); } )
                .def("replace_sequence",[] (secondary_structure::RNAStructure  & ptr, String const & seq) {
                    ptr.replace_sequence(seq); } )
                .def("get_residue",[] (secondary_structure::RNAStructure  & ptr, int num, String const & chain_id, String const & i_code) ->  secondary_structure::ResidueOP {
                    return ptr.get_residue(num, chain_id, i_code); } )
                .def("get_residue",[] (secondary_structure::RNAStructure  & ptr, util::Uuid const & uuid) ->  secondary_structure::ResidueOP {
                    return ptr.get_residue(uuid); } )
                .def("sequence",[] (secondary_structure::RNAStructure  & ptr) ->  String {
                    return ptr.sequence(); } )
                .def("dot_bracket",[] (secondary_structure::RNAStructure  & ptr) ->  String {
                    return ptr.dot_bracket(); } )
                .def("chains",[] (secondary_structure::RNAStructure  & ptr) ->  secondary_structure::ChainOPs const & {
                    return ptr.chains(); } )
                .def("residues",[] (secondary_structure::RNAStructure  & ptr) ->  secondary_structure::ResidueOPs {
                    return ptr.residues(); } )
                .def("structure",[] (secondary_structure::RNAStructure  & ptr) ->  secondary_structure::StructureOP {
                    return ptr.structure(); } )
                .def("basepairs",[] (secondary_structure::RNAStructure  & ptr) ->  secondary_structure::BasepairOPs const & {
                    return ptr.basepairs(); } )
                .def("ends",[] (secondary_structure::RNAStructure  & ptr) ->  secondary_structure::BasepairOPs const & {
                    return ptr.ends(); } )
                .def("name",[] (secondary_structure::RNAStructure  & ptr) ->  String const & {
                    return ptr.name(); } )
                .def("end_ids",[] (secondary_structure::RNAStructure  & ptr) ->  Strings const & {
                    return ptr.end_ids(); } )
                .def("name",[] (secondary_structure::RNAStructure  & ptr, String const & name) {
                    ptr.name(name); } )
                .def("path",[] (secondary_structure::RNAStructure  & ptr, String const & path) {
                    ptr.path(path); } )
                .def("end_ids",[] (secondary_structure::RNAStructure  & ptr, Strings const & end_ids) {
                    ptr.end_ids(end_ids); } )
                ;

        py::enum_<ResType>(m, "ResType")
                .value("ADE", ResType::ADE)
                .value("CYT", ResType::CYT)
                .value("GUA", ResType::GUA)
                .value("URA", ResType::URA)
                .value("NONE", ResType::NONE)

                ;

        py::class_<secondary_structure::Residue, std::shared_ptr<secondary_structure::Residue>>(m, "Residue")
                // ctors
                .def(py::init<String const &,String const &,int const &,String const &,util::Uuid const &,String const &>())
                .def(py::init<secondary_structure::Residue const &>())
                .def(py::init<String const &>())
                        // methods
                .def("to_str",[] (secondary_structure::Residue  & ptr) ->  String {
                    return ptr.to_str(); } )
                .def("name",[] (secondary_structure::Residue  & ptr) ->  String const & {
                    return ptr.name(); } )
                .def("dot_bracket",[] (secondary_structure::Residue  & ptr) ->  String const & {
                    return ptr.dot_bracket(); } )
                .def("num",[] (secondary_structure::Residue  & ptr) ->  int const & {
                    return ptr.num(); } )
                .def("chain_id",[] (secondary_structure::Residue  & ptr) ->  String const & {
                    return ptr.chain_id(); } )
                .def("i_code",[] (secondary_structure::Residue  & ptr) ->  String const & {
                    return ptr.i_code(); } )
                .def("i_code",[] (secondary_structure::Residue  & ptr, String const & code) {
                    ptr.i_code(code); } )
                .def("uuid",[] (secondary_structure::Residue  & ptr) ->  util::Uuid const & {
                    return ptr.uuid(); } )
                .def("res_type",[] (secondary_structure::Residue  & ptr) ->  secondary_structure::ResType {
                    return ptr.res_type(); } )
                .def("uuid",[] (secondary_structure::Residue  & ptr, util::Uuid const & nuuid) {
                    ptr.uuid(nuuid); } )
                .def("name",[] (secondary_structure::Residue  & ptr, String const & name) {
                    ptr.name(name); } )
                ;

        py::class_<SecondaryStructureChainGraph, std::shared_ptr<SecondaryStructureChainGraph>>(m, "SecondaryStructureChainGraph")
                // ctors
                .def(py::init<>())
                        // methods
                .def("begin",[] (SecondaryStructureChainGraph  & ptr) -> SecondaryStructureChainGraph::iterator {
                    return ptr.begin(); } )
                .def("end",[] (SecondaryStructureChainGraph  & ptr) -> SecondaryStructureChainGraph::iterator {
                    return ptr.end(); } )
                .def("begin",[] (SecondaryStructureChainGraph const & ptr) -> SecondaryStructureChainGraph::const_iterator {
                    return ptr.begin(); } )
                .def("end",[] (SecondaryStructureChainGraph const & ptr) -> SecondaryStructureChainGraph::const_iterator {
                    return ptr.end(); } )
                .def("size",[] (SecondaryStructureChainGraph  & ptr) -> size_t {
                    return ptr.size(); } )
                .def("nodes",[] (SecondaryStructureChainGraph  & ptr) -> data_structure::graph::GraphNodeOPs<NodeData> const & {
                    return ptr.nodes(); } )
                .def("add_chain",[] (SecondaryStructureChainGraph  & ptr, NodeData const & data, int parent_index, int orphan) -> int {
                    return ptr.add_chain(data, parent_index, orphan); } )
                .def("get_node_by_res",[] (SecondaryStructureChainGraph  & ptr, secondary_structure::ResidueOP const & res) -> int {
                    return ptr.get_node_by_res(res); } )
                .def("pair_res",[] (SecondaryStructureChainGraph  & ptr, int n_i, int n_j) {
                    ptr.pair_res(n_i, n_j); } )
                ;

        py::class_<SecondaryStructureTree, std::shared_ptr<SecondaryStructureTree>>(m, "SecondaryStructureTree")
                // ctors
                .def(py::init<>())
                        // methods
                .def("begin",[] (SecondaryStructureTree  & ptr) -> SecondaryStructureTree::iterator {
                    return ptr.begin(); } )
                .def("end",[] (SecondaryStructureTree  & ptr) -> SecondaryStructureTree::iterator {
                    return ptr.end(); } )
                .def("begin",[] (SecondaryStructureTree const & ptr) -> SecondaryStructureTree::const_iterator {
                    return ptr.begin(); } )
                .def("end",[] (SecondaryStructureTree const & ptr) -> SecondaryStructureTree::const_iterator {
                    return ptr.end(); } )
                .def("size",[] (SecondaryStructureTree  & ptr) -> size_t {
                    return ptr.size(); } )
                .def("add_motif",[] (SecondaryStructureTree  & ptr, secondary_structure::MotifOP const & m, int parent_index, int parent_end_index) -> int {
                    return ptr.add_motif(m, parent_index, parent_end_index); } )
                ;
/*
        py::class_<SequenceConstraint, std::shared_ptr<SequenceConstraint>>(m, "SequenceConstraint")
		// ctors
		.def(py::init<>())
		// methods
		.def("clone",[] (SequenceConstraint const & ptr) -> SequenceConstraint * {
		 return ptr.clone(); } )
		.def("violates_constraint",[] (SequenceConstraint  & ptr, PoseOP p) -> bool {
		 return ptr.violates_constraint(p); } )
		.def("violations",[] (SequenceConstraint  & ptr, PoseOP ) -> int {
		 return ptr.violations(); } )
		;
*/
        py::class_<SequenceConstraints, std::shared_ptr<SequenceConstraints>>(m, "SequenceConstraints")
                // ctors
                .def(py::init<>())
                        // methods
                .def("add_sequence_constraint",[] (SequenceConstraints  & ptr, secondary_structure::SequenceConstraintOP seq_constraint) {
                    ptr.add_sequence_constraint(seq_constraint); } )
                .def("add_disallowed_sequence",[] (SequenceConstraints  & ptr, String const & seq) {
                    ptr.add_disallowed_sequence(seq); } )
                .def("add_gc_helix_stretch_limit",[] (SequenceConstraints  & ptr, int length) {
                    ptr.add_gc_helix_stretch_limit(length); } )
                .def("violations",[] (SequenceConstraints  & ptr, secondary_structure::PoseOP p) -> Ints const & {
                    return ptr.violations(p); } )
                .def("num_constraints",[] (SequenceConstraints  & ptr) -> size_t {
                    return ptr.num_constraints(); } )
                ;

        py::class_<secondary_structure::Structure, std::shared_ptr<secondary_structure::Structure>>(m, "SecondaryStructureStructure")
                // ctors
                .def(py::init<secondary_structure::ChainOPs const &>())
                .def(py::init<String const &,String const &>())
                .def(py::init<secondary_structure::Structure const &>())
                .def(py::init<String const &>())
                        // methods
                .def("residues",[] (secondary_structure::Structure  & ptr) ->  secondary_structure::ResidueOPs {
                    return ptr.residues(); } )
                .def("sequence",[] (secondary_structure::Structure  & ptr) ->  String {
                    return ptr.sequence(); } )
                .def("dot_bracket",[] (secondary_structure::Structure  & ptr) ->  String {
                    return ptr.dot_bracket(); } )
                .def("get_residue",[] (secondary_structure::Structure  & ptr, int const & num, String const & chain_id, String const & i_code) -> secondary_structure::ResidueOP  {
                    return ptr.get_residue(num, chain_id, i_code); } )
                .def("get_residue",[] (secondary_structure::Structure  & ptr, int const & num, String const & chain_id, String const & i_code) ->secondary_structure:: ResidueOP  {
                    return ptr.get_residue(num, chain_id, i_code); } )
                .def("to_str",[] (secondary_structure::Structure  & ptr) -> String {
                    return ptr.to_str(); } )
                .def("chains",[] (secondary_structure::Structure  & ptr) ->  secondary_structure::ChainOPs const & {
                    return ptr.chains(); } )
                ;

/*
        py::class_<TreePreNode, std::shared_ptr<TreePreNode>>(m, "TreePreNode")
		// ctors
		.def(py::init<>())
		.def(py::init<MotifOP const &,int,int>())
		// public attributes
		.def_readwrite("m", &TreePreNode::m)
		.def_readwrite("parent_index", &TreePreNode::parent_index)
		.def_readwrite("parent_end_index", &TreePreNode::parent_end_index)
		;
        py::class_<res_less_than_key, std::shared_ptr<res_less_than_key>>(m, "res_less_than_key")
		// operators
		.def(py::self () secondary_structure::ResidueOP)
		;
*/
        // enums
        py::enum_<NodeType>(m,"NodeType")
                .value("UNPAIRED",NodeType::UNPAIRED)
                .value("PAIRED",NodeType::PAIRED)
                ;
    }
}

#endif // PYBIND11_SECONDARY_STRUCTURE_HPP