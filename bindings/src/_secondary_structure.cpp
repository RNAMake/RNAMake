// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/cast.h>
#include <pybind11/operators.h>
#include <memory>

// base includes
#include <base/application.hpp>
#include <base/backtrace.h>
#include <base/cl_option.h>
#include <base/command_line_parser.hpp>
#include <base/env_manager.h>
#include <base/exception.h>
#include <base/file_io.h>
#include <base/log.h>
#include <base/option.h>
#include <base/settings.h>
#include <base/string.h>
#include <base/sys_interface.h>
#include <base/vector_container.h>

// data_structure includes
#include <data_structure/graph.h>
#include <data_structure/graph_adjacency_list.h>
#include <data_structure/graph_base.h>
#include <data_structure/graph_iter_list.h>

// data_structure::graph includes
#include <data_structure/graph/graph.h>
#include <data_structure/graph/graph_node.fwd.h>
#include <data_structure/graph/graph_node.h>

// data_structure::tree includes
#include <data_structure/tree/tree.h>
#include <data_structure/tree/tree_node.h>

// eternabot includes
#include <eternabot/feature_generator.h>
#include <eternabot/scorer.h>
#include <eternabot/sequence_designer.h>
#include <eternabot/strategy.h>
#include <eternabot/strategy/a_basic_test.h>
#include <eternabot/strategy/berex_test.h>
#include <eternabot/strategy/clear_plot.h>
#include <eternabot/strategy/direction_of_gc.h>
#include <eternabot/strategy/num_of_yellow.h>

// io includes
#include <util/csv.h>

// io::detail includes
#include <util/csv.h>

// io::error includes
#include <util/csv.h>

// math includes
#include <math/euler.h>
#include <math/hashing.h>
#include <math/numerical.h>
#include <math/quaternion.h>
#include <math/stats.h>
#include <math/transform.h>
#include <math/xyz_matrix.h>
#include <math/xyz_vector.h>

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

// motif_data_structure includes
#include <motif_data_structure/motif_connection.h>
#include <motif_data_structure/motif_graph.h>
#include <motif_data_structure/motif_merger.h>
#include <motif_data_structure/motif_state_ensemble_graph.h>
#include <motif_data_structure/motif_state_ensemble_tree.h>
#include <motif_data_structure/motif_state_graph.hpp>
#include <motif_data_structure/motif_state_node.hpp>
#include <motif_data_structure/motif_state_tree.h>
#include <motif_data_structure/motif_topology.h>
#include <motif_data_structure/motif_tree.h>

// motif_search includes
#include <motif_search/motif_state_monte_carlo.h>
#include <motif_search/problem.h>
#include <motif_search/search.h>
#include <motif_search/solution_filter.h>
#include <motif_search/solution_topology.h>

// motif_search::exhaustive includes
#include <motif_search/exhaustive/motif_state_enumerator.h>
#include <motif_search/exhaustive/scorer.h>
#include <motif_search/exhaustive/search.h>

// motif_search::monte_carlo includes
#include <motif_search/monte_carlo/scorer.h>
#include <motif_search/monte_carlo/search.h>

// motif_search::path_finding includes
#include <motif_search/path_finding/node.h>
#include <motif_search/path_finding/scorer.h>
#include <motif_search/path_finding/search.h>
#include <motif_search/path_finding/selector.h>

// motif_tools includes
#include <motif_tools/segmenter.h>

// resources includes
#include <resources/added_motif_library.h>
#include <resources/motif_ensemble_sqlite_connection.h>
#include <resources/motif_sqlite_connection.h>
#include <resources/motif_sqlite_library.h>
#include <resources/motif_state_ensemble_sqlite_library.h>
#include <resources/motif_state_sqlite_library.h>
#include <resources/resource_manager.h>
#include <resources/sqlite_library.h>

// secondary_structure includes
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

// sequence_optimization includes
#include <sequence_optimization/sequence_optimizer_3d.hpp>

// structure includes
#include <structure/atom.h>
#include <structure/basepair.h>
#include <structure/basepair_state.h>
#include <structure/beads.h>
#include <structure/chain.h>
#include <structure/cif_parser.h>
#include <structure/close_chain.h>
#include <structure/is_equal.h>
#include <structure/pdb_parser.h>
#include <structure/residue.h>
#include <structure/residue_type.h>
#include <structure/residue_type_set.h>
#include <structure/residue_type_set_manager.h>
#include <structure/rna_structure.h>
#include <structure/structure.h>

// thermo_fluctuation includes
#include <thermo_fluctuation/thermo_fluc_relax.h>
#include <thermo_fluctuation/thermo_fluc_sampler.h>
#include <thermo_fluctuation/thermo_fluc_scorer.h>
#include <thermo_fluctuation/thermo_fluc_simulation.h>

// thermo_fluctuation::graph includes
#include <thermo_fluctuation/graph/sampler.h>
#include <thermo_fluctuation/graph/scorer.h>
#include <thermo_fluctuation/graph/simulation.h>

// thermo_fluctuation::graph::logging includes
#include <thermo_fluctuation/graph/logging.h>

// thermo_fluctuation::graph::sterics includes
#include <thermo_fluctuation/graph/sterics.h>

// util includes
#include <util/basic_io.hpp>
#include <util/cartesian_product.h>
#include <util/monte_carlo.h>
#include <util/motif_type.h>
#include <util/random_number_generator.h>
#include <util/sqlite3_connection.h>
#include <util/steric_lookup.hpp>
#include <util/uuid.h>
#include <util/x3dna.h>
using namespace secondary_structure;

// Bindings start

namespace py = pybind11;

PYBIND11_MODULE(_secondary_structure,m) {

	// free functions

    m.def("assign_end_id", [] (RNAStructureOP const & ss, BasepairOP const & end) -> String {
		 return assign_end_id(ss, end); },
	py::arg("ss"),
	py::arg("end") 
    );

    m.def("convert_res_name_to_type", [] (char c) -> ResType {
		 return convert_res_name_to_type(c); },
	py::arg("c") 
    );

    m.def("fill_basepairs_in_ss", [] (PoseOP & ss) {
		fill_basepairs_in_ss(ss); },
	py::arg("ss") 
    );

    m.def("find_gc_helix_stretches", [] (PoseOP p, int length) -> int {
		 return find_gc_helix_stretches(p, length); },
	py::arg("p"),
	py::arg("length") 
    );

    m.def("find_longest_gc_helix_stretch", [] (PoseOP p) -> int {
		 return find_longest_gc_helix_stretch(p); },
	py::arg("p") 
    );

    m.def("find_res_types_in_pose", [] (PoseOP p, ResTypes const & residue_types) -> int {
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

    m.def("tree_from_pose", [] (PoseOP const & p) -> SecondaryStructureTreeOP {
		 return tree_from_pose(p); },
	py::arg("p") 
    );
	// classes

        py::class_<Basepair, std::shared_ptr<Basepair>>(m, "Basepair")
		// ctors
		.def(py::init<ResidueOP const &,ResidueOP const &,util::Uuid const &>())
		// methods
		.def("name",[] (Basepair  & ptr) ->  String {
		 return ptr.name(); } )
		.def("partner",[] (Basepair  & ptr, secondary_structure::ResidueOP const & r) ->  ResidueOP {
		 return ptr.partner(r); } )
		.def("res1",[] (Basepair  & ptr) ->  ResidueOP & {
		 return ptr.res1(); } )
		.def("res2",[] (Basepair  & ptr) ->  ResidueOP & {
		 return ptr.res2(); } )
		.def("res1",[] (Basepair const & ptr) ->  ResidueOP const & {
		 return ptr.res1(); } )
		.def("res2",[] (Basepair const & ptr) ->  ResidueOP const & {
		 return ptr.res2(); } )
		.def("uuid",[] (Basepair  & ptr) ->  util::Uuid const & {
		 return ptr.uuid(); } )
		;

        py::class_<Chain, std::shared_ptr<Chain>>(m, "Chain")
		// ctors
		.def(py::init<>())
		.def(py::init<ResidueOPs const &>())
		.def(py::init<Chain const &>())
		.def(py::init<String const &>())
		// methods
		.def("begin",[] (Chain  & ptr) -> ResidueOPs::iterator {
		 return ptr.begin(); } )
		.def("end",[] (Chain  & ptr) -> ResidueOPs::iterator {
		 return ptr.end(); } )
		.def("begin",[] (Chain const & ptr) -> ResidueOPs::const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (Chain const & ptr) -> ResidueOPs::const_iterator {
		 return ptr.end(); } )
		.def("first",[] (Chain  & ptr) ->  ResidueOP const & {
		 return ptr.first(); } )
		.def("last",[] (Chain  & ptr) ->  ResidueOP const & {
		 return ptr.last(); } )
		.def("sequence",[] (Chain  & ptr) ->  String {
		 return ptr.sequence(); } )
		.def("dot_bracket",[] (Chain  & ptr) ->  String {
		 return ptr.dot_bracket(); } )
		.def("to_str",[] (Chain  & ptr) ->  String {
		 return ptr.to_str(); } )
		.def("length",[] (Chain  & ptr) ->  int {
		 return ptr.length(); } )
		.def("residues",[] (Chain  & ptr) ->  ResidueOPs const & {
		 return ptr.residues(); } )
		;

        py::class_<DisallowedSequence, std::shared_ptr<DisallowedSequence>>(m, "DisallowedSequence")
		// ctors
		.def(py::init<String const &>())
		// methods
		.def("clone",[] (DisallowedSequence const & ptr) -> SequenceConstraint * {
		 return ptr.clone(); } )
		.def("violations",[] (DisallowedSequence  & ptr, PoseOP p) -> int {
		 return ptr.violations(p); } )
		// inherited methods
		.def("clone",[] (SequenceConstraint const & ptr) -> SequenceConstraint * {
		 return ptr.clone(); } )
		.def("violates_constraint",[] (SequenceConstraint  & ptr, PoseOP p) -> bool {
		 return ptr.violates_constraint(p); } )
		.def("violations",[] (SequenceConstraint  & ptr, PoseOP p) -> int {
		 return ptr.violations(p); } )
		;

        py::class_<GCHelixStretchLimit, std::shared_ptr<GCHelixStretchLimit>>(m, "GCHelixStretchLimit")
		// ctors
		.def(py::init<int>())
		// methods
		.def("clone",[] (GCHelixStretchLimit const & ptr) -> SequenceConstraint * {
		 return ptr.clone(); } )
		.def("violations",[] (GCHelixStretchLimit  & ptr, PoseOP p) -> int {
		 return ptr.violations(p); } )
		// inherited methods
		.def("clone",[] (SequenceConstraint const & ptr) -> SequenceConstraint * {
		 return ptr.clone(); } )
		.def("violates_constraint",[] (SequenceConstraint  & ptr, PoseOP p) -> bool {
		 return ptr.violates_constraint(p); } )
		.def("violations",[] (SequenceConstraint  & ptr, PoseOP p) -> int {
		 return ptr.violations(p); } )
		;

        py::class_<Motif, std::shared_ptr<Motif>>(m, "Motif")
		// ctors
		.def(py::init<>())
		.def(py::init<StructureOP const &,BasepairOPs const &,BasepairOPs const &>())
		.def(py::init<StructureOP const &,BasepairOPs const &,BasepairOPs const &,Strings const &,String const &,String const &,float>())
		.def(py::init<Motif const &>())
		.def(py::init<String const &>())
		// methods
		.def("to_str",[] (Motif  & ptr) -> String {
		 return ptr.to_str(); } )
		.def("mtype",[] (Motif  & ptr) ->  util::MotifType const & {
		 return ptr.mtype(); } )
		.def("id",[] (Motif  & ptr) ->  util::Uuid const & {
		 return ptr.id(); } )
		.def("mtype",[] (Motif  & ptr, util::MotifType const & mtype) {
		ptr.mtype(mtype); } )
		.def("id",[] (Motif  & ptr, util::Uuid const & uuid) {
		ptr.id(uuid); } )
		// inherited methods
		.def("get_basepair",[] (RNAStructure  & ptr, util::Uuid const & bp_uuid) -> BasepairOPs {
		 return ptr.get_basepair(bp_uuid); } )
		.def("get_basepair",[] (RNAStructure  & ptr, util::Uuid const & bp_uuid) -> BasepairOPs {
		 return ptr.get_basepair(bp_uuid); } )
		.def("get_basepair",[] (RNAStructure  & ptr, util::Uuid const & bp_uuid) -> BasepairOPs {
		 return ptr.get_basepair(bp_uuid); } )
		.def("get_basepair",[] (RNAStructure  & ptr, util::Uuid const & bp_uuid) -> BasepairOPs {
		 return ptr.get_basepair(bp_uuid); } )
		.def("get_end",[] (RNAStructure  & ptr, String const & name) -> BasepairOP {
		 return ptr.get_end(name); } )
		.def("replace_sequence",[] (RNAStructure  & ptr, String const & seq) {
		ptr.replace_sequence(seq); } )
		.def("get_residue",[] (RNAStructure  & ptr, int num, String const & chain_id, String const & i_code) ->  ResidueOP {
		 return ptr.get_residue(num, chain_id, i_code); } )
		.def("get_residue",[] (RNAStructure  & ptr, util::Uuid const & uuid) ->  ResidueOP {
		 return ptr.get_residue(uuid); } )
		.def("sequence",[] (RNAStructure  & ptr) ->  String {
		 return ptr.sequence(); } )
		.def("dot_bracket",[] (RNAStructure  & ptr) ->  String {
		 return ptr.dot_bracket(); } )
		.def("chains",[] (RNAStructure  & ptr) ->  ChainOPs const & {
		 return ptr.chains(); } )
		.def("residues",[] (RNAStructure  & ptr) ->  ResidueOPs {
		 return ptr.residues(); } )
		.def("structure",[] (RNAStructure  & ptr) ->  StructureOP {
		 return ptr.structure(); } )
		.def("basepairs",[] (RNAStructure  & ptr) ->  BasepairOPs const & {
		 return ptr.basepairs(); } )
		.def("ends",[] (RNAStructure  & ptr) ->  BasepairOPs const & {
		 return ptr.ends(); } )
		.def("name",[] (RNAStructure  & ptr) ->  String const & {
		 return ptr.name(); } )
		.def("end_ids",[] (RNAStructure  & ptr) ->  Strings const & {
		 return ptr.end_ids(); } )
		.def("name",[] (RNAStructure  & ptr, String const & name) {
		ptr.name(name); } )
		.def("path",[] (RNAStructure  & ptr, String const & path) {
		ptr.path(path); } )
		.def("end_ids",[] (RNAStructure  & ptr, Strings const & end_ids) {
		ptr.end_ids(end_ids); } )
		;

        py::class_<NodeData, std::shared_ptr<NodeData>>(m, "NodeData")
		// ctors
		.def(py::init<>())
		.def(py::init<ResidueOPs const &,NodeType const &>())
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
		.def("parse_to_motif",[] (Parser  & ptr, String const & sequence, String const & dot_bracket) -> MotifOP  {
		 return ptr.parse_to_motif(sequence, dot_bracket); } )
		.def("parse_to_pose",[] (Parser  & ptr, String const & sequence, String const & dot_bracket) -> PoseOP {
		 return ptr.parse_to_pose(sequence, dot_bracket); } )
		.def("reset",[] (Parser  & ptr) {
		ptr.reset(); } )
		;

        py::class_<Pose, std::shared_ptr<Pose>>(m, "Pose")
		// ctors
		.def(py::init<>())
		.def(py::init<StructureOP const &,BasepairOPs const &,BasepairOPs const &>())
		.def(py::init<StructureOP const &,BasepairOPs const &,BasepairOPs const &,MotifOPs const &>())
		.def(py::init<RNAStructureOP const &,MotifOPs const &>())
		// methods
		.def("helices",[] (Pose  & ptr) -> MotifOPs const & {
		 return ptr.helices(); } )
		.def("motifs",[] (Pose const & ptr) -> MotifOPs const & {
		 return ptr.motifs(); } )
		.def("motif",[] (Pose  & ptr, util::Uuid const & uuid) -> MotifOP {
		 return ptr.motif(uuid); } )
		.def("replace_sequence",[] (Pose  & ptr, String const & seq) {
		ptr.replace_sequence(seq); } )
		.def("update_motif",[] (Pose  & ptr, util::Uuid const & uuid) {
		ptr.update_motif(uuid); } )
		// inherited methods
		.def("get_basepair",[] (RNAStructure  & ptr, util::Uuid const & bp_uuid) -> BasepairOPs {
		 return ptr.get_basepair(bp_uuid); } )
		.def("get_basepair",[] (RNAStructure  & ptr, util::Uuid const & bp_uuid) -> BasepairOPs {
		 return ptr.get_basepair(bp_uuid); } )
		.def("get_basepair",[] (RNAStructure  & ptr, util::Uuid const & bp_uuid) -> BasepairOPs  {
		 return ptr.get_basepair(bp_uuid); } )
		.def("get_basepair",[] (RNAStructure  & ptr, util::Uuid const & bp_uuid) -> BasepairOPs {
		 return ptr.get_basepair(bp_uuid); } )
		.def("get_end",[] (RNAStructure  & ptr, String const & name) -> BasepairOP {
		 return ptr.get_end(name); } )
		.def("replace_sequence",[] (RNAStructure  & ptr, String const & seq) {
		ptr.replace_sequence(seq); } )
		.def("get_residue",[] (RNAStructure  & ptr, int num, String const & chain_id, String const & i_code) ->  ResidueOP {
		 return ptr.get_residue(num, chain_id, i_code); } )
		.def("get_residue",[] (RNAStructure  & ptr, util::Uuid const & uuid) ->  ResidueOP {
		 return ptr.get_residue(uuid); } )
		.def("sequence",[] (RNAStructure  & ptr) ->  String {
		 return ptr.sequence(); } )
		.def("dot_bracket",[] (RNAStructure  & ptr) ->  String {
		 return ptr.dot_bracket(); } )
		.def("chains",[] (RNAStructure  & ptr) ->  ChainOPs const & {
		 return ptr.chains(); } )
		.def("residues",[] (RNAStructure  & ptr) ->  ResidueOPs {
		 return ptr.residues(); } )
		.def("structure",[] (RNAStructure  & ptr) ->  StructureOP {
		 return ptr.structure(); } )
		.def("basepairs",[] (RNAStructure  & ptr) ->  BasepairOPs const & {
		 return ptr.basepairs(); } )
		.def("ends",[] (RNAStructure  & ptr) ->  BasepairOPs const & {
		 return ptr.ends(); } )
		.def("name",[] (RNAStructure  & ptr) ->  String const & {
		 return ptr.name(); } )
		.def("end_ids",[] (RNAStructure  & ptr) ->  Strings const & {
		 return ptr.end_ids(); } )
		.def("name",[] (RNAStructure  & ptr, String const & name) {
		ptr.name(name); } )
		.def("path",[] (RNAStructure  & ptr, String const & path) {
		ptr.path(path); } )
		.def("end_ids",[] (RNAStructure  & ptr, Strings const & end_ids) {
		ptr.end_ids(end_ids); } )
		;

        py::class_<RNAStructure, std::shared_ptr<RNAStructure>>(m, "RNAStructure")
		// ctors
		.def(py::init<>())
		.def(py::init<StructureOP const &,BasepairOPs const &,BasepairOPs const &>())
		.def(py::init<StructureOP const &,BasepairOPs const &,BasepairOPs const &,Strings const &,String const &,String const &,float>())
		.def(py::init<RNAStructure const &>())
		// methods
		.def("get_basepair",[] (RNAStructure  & ptr, util::Uuid const & bp_uuid) -> BasepairOPs {
		 return ptr.get_basepair(bp_uuid); } )
		.def("get_basepair",[] (RNAStructure  & ptr, util::Uuid const & bp_uuid) -> BasepairOPs {
		 return ptr.get_basepair(bp_uuid); } )
		.def("get_basepair",[] (RNAStructure  & ptr, util::Uuid const & bp_uuid) -> BasepairOPs {
		 return ptr.get_basepair(bp_uuid); } )
		.def("get_basepair",[] (RNAStructure  & ptr, util::Uuid const & bp_uuid) -> BasepairOPs {
		 return ptr.get_basepair(bp_uuid); } )
		.def("get_end",[] (RNAStructure  & ptr, String const & name) -> BasepairOP  {
		 return ptr.get_end(name); } )
		.def("replace_sequence",[] (RNAStructure  & ptr, String const & seq) {
		ptr.replace_sequence(seq); } )
		.def("get_residue",[] (RNAStructure  & ptr, int num, String const & chain_id, String const & i_code) ->  ResidueOP {
		 return ptr.get_residue(num, chain_id, i_code); } )
		.def("get_residue",[] (RNAStructure  & ptr, util::Uuid const & uuid) ->  ResidueOP {
		 return ptr.get_residue(uuid); } )
		.def("sequence",[] (RNAStructure  & ptr) ->  String {
		 return ptr.sequence(); } )
		.def("dot_bracket",[] (RNAStructure  & ptr) ->  String {
		 return ptr.dot_bracket(); } )
		.def("chains",[] (RNAStructure  & ptr) ->  ChainOPs const & {
		 return ptr.chains(); } )
		.def("residues",[] (RNAStructure  & ptr) ->  ResidueOPs {
		 return ptr.residues(); } )
		.def("structure",[] (RNAStructure  & ptr) ->  StructureOP {
		 return ptr.structure(); } )
		.def("basepairs",[] (RNAStructure  & ptr) ->  BasepairOPs const & {
		 return ptr.basepairs(); } )
		.def("ends",[] (RNAStructure  & ptr) ->  BasepairOPs const & {
		 return ptr.ends(); } )
		.def("name",[] (RNAStructure  & ptr) ->  String const & {
		 return ptr.name(); } )
		.def("end_ids",[] (RNAStructure  & ptr) ->  Strings const & {
		 return ptr.end_ids(); } )
		.def("name",[] (RNAStructure  & ptr, String const & name) {
		ptr.name(name); } )
		.def("path",[] (RNAStructure  & ptr, String const & path) {
		ptr.path(path); } )
		.def("end_ids",[] (RNAStructure  & ptr, Strings const & end_ids) {
		ptr.end_ids(end_ids); } )
		;

        py::class_<ResType, std::shared_ptr<ResType>>(m, "ResType")
		;

        py::class_<Residue, std::shared_ptr<Residue>>(m, "Residue")
		// ctors
		.def(py::init<String const &,String const &,int const &,String const &,util::Uuid const &,String const &>())
		.def(py::init<Residue const &>())
		.def(py::init<String const &>())
		// methods
		.def("to_str",[] (Residue  & ptr) ->  String {
		 return ptr.to_str(); } )
		.def("name",[] (Residue  & ptr) ->  String const & {
		 return ptr.name(); } )
		.def("dot_bracket",[] (Residue  & ptr) ->  String const & {
		 return ptr.dot_bracket(); } )
		.def("num",[] (Residue  & ptr) ->  int const & {
		 return ptr.num(); } )
		.def("chain_id",[] (Residue  & ptr) ->  String const & {
		 return ptr.chain_id(); } )
		.def("i_code",[] (Residue  & ptr) ->  String const & {
		 return ptr.i_code(); } )
		.def("i_code",[] (Residue  & ptr, String const & code) {
		ptr.i_code(code); } )
		.def("uuid",[] (Residue  & ptr) ->  util::Uuid const & {
		 return ptr.uuid(); } )
		.def("res_type",[] (Residue  & ptr) ->  ResType {
		 return ptr.res_type(); } )
		.def("uuid",[] (Residue  & ptr, util::Uuid const & nuuid) {
		ptr.uuid(nuuid); } )
		.def("name",[] (Residue  & ptr, String const & name) {
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
		.def("violations",[] (SequenceConstraints  & ptr, PoseOP p) -> Ints const & {
		 return ptr.violations(p); } )
		.def("num_constraints",[] (SequenceConstraints  & ptr) -> size_t {
		 return ptr.num_constraints(); } )
		;

        py::class_<Structure, std::shared_ptr<Structure>>(m, "Structure")
		// ctors
		.def(py::init<ChainOPs const &>())
		.def(py::init<String const &,String const &>())
		.def(py::init<Structure const &>())
		.def(py::init<String const &>())
		// methods
		.def("residues",[] (Structure  & ptr) ->  ResidueOPs {
		 return ptr.residues(); } )
		.def("sequence",[] (Structure  & ptr) ->  String {
		 return ptr.sequence(); } )
		.def("dot_bracket",[] (Structure  & ptr) ->  String {
		 return ptr.dot_bracket(); } )
		.def("get_residue",[] (Structure  & ptr, int const & num, String const & chain_id, String const & i_code) -> ResidueOP  {
		 return ptr.get_residue(num, chain_id, i_code); } )
		.def("get_residue",[] (Structure  & ptr, int const & num, String const & chain_id, String const & i_code) -> ResidueOP  {
		 return ptr.get_residue(num, chain_id, i_code); } )
		.def("to_str",[] (Structure  & ptr) -> String {
		 return ptr.to_str(); } )
		.def("chains",[] (Structure  & ptr) ->  ChainOPs const & {
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