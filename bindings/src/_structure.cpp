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
using namespace structure;

// Bindings start

namespace py = pybind11;

PYBIND11_MODULE(_structure,m) {

	// free functions

    m.def("are_atom_vectors_equal", [] (AtomOPs const & atoms_1, AtomOPs const & atoms_2, float tol) -> bool {
		 return are_atom_vectors_equal(atoms_1, atoms_2, tol); },
	py::arg("atoms_1"),
	py::arg("atoms_2"),
	py::arg("tol") 
    );

    m.def("are_atoms_equal", [] (AtomOP const & a1, AtomOP const & a2, float tol) -> bool {
		 return are_atoms_equal(a1, a2, tol); },
	py::arg("a1"),
	py::arg("a2"),
	py::arg("tol") 
    );

    m.def("are_basepair_states_equal", [] (BasepairState const & a, BasepairState const & b) -> int {
		 return are_basepair_states_equal(a, b); },
	py::arg("a"),
	py::arg("b") 
    );

    m.def("are_chains_equal", [] (ChainOP const & c1, ChainOP const & c2, int check_uuids) -> bool {
		 return are_chains_equal(c1, c2, check_uuids); },
	py::arg("c1"),
	py::arg("c2"),
	py::arg("check_uuids") = 1 
    );

    m.def("are_residues_equal", [] (ResidueOP const & r1, ResidueOP const & r2, int check_uuids) -> bool {
		 return are_residues_equal(r1, r2, check_uuids); },
	py::arg("r1"),
	py::arg("r2"),
	py::arg("check_uuids") = 1 
    );

    m.def("are_structures_equal", [] (StructureOP const & s1, StructureOP const & s2, int check_uuids) -> bool {
		 return are_structures_equal(s1, s2, check_uuids); },
	py::arg("s1"),
	py::arg("s2"),
	py::arg("check_uuids") = 1 
    );

    m.def("axis_angle_to_rot_matrix", [] (float angle, math::Vector const & al) -> math::Matrix {
		 return axis_angle_to_rot_matrix(angle, al); },
	py::arg("angle"),
	py::arg("al") 
    );

    m.def("center", [] (AtomOPs const & atoms) -> math::Point {
		 return center(atoms); },
	py::arg("atoms") 
    );

    m.def("close_chain", [] (ChainOP chain) {
		close_chain(chain); },
	py::arg("chain") 
    );

    m.def("close_torsion", [] (int which_dir, AtomOPs const & parent_atoms, AtomOPs daughter_atoms, AtomOPs const & match_atoms_1, AtomOPs const & match_atoms_2) {
		close_torsion(which_dir, parent_atoms, daughter_atoms, match_atoms_1, match_atoms_2); },
	py::arg("which_dir"),
	py::arg("parent_atoms"),
	py::arg("daughter_atoms"),
	py::arg("match_atoms_1"),
	py::arg("match_atoms_2") 
    );

    m.def("connect_residues_into_chains", [] (ResidueOPs & residues, ChainOPs & chains) {
		connect_residues_into_chains(residues, chains); },
	py::arg("residues"),
	py::arg("chains") 
    );

    m.def("create_coord_system", [] (AtomOPs const & atoms) -> math::Matrix {
		 return create_coord_system(atoms); },
	py::arg("atoms") 
    );

    m.def("end_from_basepairs", [] (StructureOP const & s, BasepairOPs const & bps) {//-> std::shared_ptr<BasepairOPs> {
		 return *end_from_basepairs(s, bps); },
	py::arg("s"),
	py::arg("bps") 
    );

    m.def("frame_distance", [] (BasepairStateOP const & current, BasepairStateOP const & end, BasepairStateOP const & endflip) ->  const float {
		 return frame_distance(current, end, endflip); },
	py::arg("current"),
	py::arg("end"),
	py::arg("endflip") 
    );

    m.def("get_bpstate_rotation_diff", [] (BasepairState const & bp1, BasepairState const & bp2) -> float {
		 return get_bpstate_rotation_diff(bp1, bp2); },
	py::arg("bp1"),
	py::arg("bp2") 
    );

    m.def("get_projection", [] (math::Point const & coord, math::Point const & current_pos, math::Vector const & projection_axis) -> math::Vector {
		 return get_projection(coord, current_pos, projection_axis); },
	py::arg("coord"),
	py::arg("current_pos"),
	py::arg("projection_axis") 
    );

    m.def("get_ref_bp_state", [] () -> BasepairState {
		 return get_ref_bp_state(); } 
    );

    m.def("get_res_ref_frame", [] (ResidueOP r) -> math::Matrix {
		 return get_res_ref_frame(r); },
	py::arg("r") 
    );

    m.def("gu_bp", [] (BasepairOP const & bp) -> bool {
		 return gu_bp(bp); },
	py::arg("bp") 
    );

    m.def("m_dot_v", [] (math::Matrix const & m, math::Vector const & v) -> math::Vector {
		 return m_dot_v(m, v); },
	py::arg("m"),
	py::arg("v") 
    );

    m.def("new_score_function", [] (BasepairStateOP const & current, BasepairStateOP const & end, BasepairStateOP const & endflip) ->  const float {
		 return new_score_function(current, end, endflip); },
	py::arg("current"),
	py::arg("end"),
	py::arg("endflip") 
    );

    m.def("new_score_function_new", [] (BasepairStateOP const & current, BasepairStateOP const & end, BasepairStateOP const & endflip) ->  const float {
		 return new_score_function_new(current, end, endflip); },
	py::arg("current"),
	py::arg("end"),
	py::arg("endflip") 
    );

    m.def("replace_missing_phosphate_backbone", [] (ResidueOP r, ResidueOP r_template) {
		replace_missing_phosphate_backbone(r, r_template); },
	py::arg("r"),
	py::arg("r_template") 
    );

    m.def("str_to_basepairstate", [] (String const & s) -> BasepairState {
		 return str_to_basepairstate(s); },
	py::arg("s") 
    );

    m.def("subselect_basepairs_with_res", [] (ResidueOPs const & res, BasepairOPs const & all_bps) {
		 return *subselect_basepairs_with_res(res, all_bps); },
	py::arg("res"),
	py::arg("all_bps") 
    );

    m.def("to_radians", [] (float degrees) -> float {
		 return to_radians(degrees); },
	py::arg("degrees") 
    );

    m.def("virtual_atom", [] (String const & name, float l, float theta, float phi, AtomOPs const & parent_atoms) -> AtomOP {
		 return virtual_atom(name, l, theta, phi, parent_atoms); },
	py::arg("name"),
	py::arg("l"),
	py::arg("theta"),
	py::arg("phi"),
	py::arg("parent_atoms") 
    );

    m.def("wc_bp", [] (BasepairOP const & bp) -> bool {
		 return wc_bp(bp); },
	py::arg("bp") 
    );
	// classes

        py::class_<Atom, std::shared_ptr<Atom>>(m, "Atom")
		// ctors
		.def(py::init<String const &,math::Point const &>())
		.def(py::init<String const &>())
		.def(py::init<Atom const &>())
		// methods
		.def("to_str",[] (Atom  & ptr) -> String {
		 return ptr.to_str(); } )
		.def("to_pdb_str",[] (Atom  & ptr, int num) -> String {
		 return ptr.to_pdb_str(num); } )
		.def("move",[] (Atom  & ptr, math::Point const & p) {
		ptr.move(p); } )
		.def("transform",[] (Atom  & ptr, math::Transform const & t) {
		ptr.transform(t); } )
		.def("name",[] (Atom const & ptr) ->  String const & {
		 return ptr.name(); } )
		.def("coords",[] (Atom const & ptr) ->  math::Point const {
		 return ptr.coords(); } )
		.def("coords",[] (Atom  & ptr, math::Point const & ncoords) {
		ptr.coords(ncoords); } )
		.def("name",[] (Atom  & ptr, String const & nname) {
		ptr.name(nname); } )
		;

        py::class_<Basepair, std::shared_ptr<Basepair>>(m, "Basepair")
		// ctors
		.def(py::init<>())
		.def(py::init<ResidueOP &,ResidueOP &,math::Matrix const &,String const &>())
		// methods
		.def("copy",[] (Basepair  & ptr) -> Basepair {
		 return ptr.copy(); } )
		.def("move",[] (Basepair  & ptr, math::Point const & p) {
		ptr.move(p); } )
		.def("transform",[] (Basepair  & ptr, math::Matrix const & r, math::Vector const & t, math::Point & dummy) {
		ptr.transform(r, t, dummy); } )
		.def("transform",[] (Basepair  & ptr, math::Matrix const & r, math::Vector const & t) {
		ptr.transform(r, t); } )
		.def("state",[] (Basepair  & ptr) ->  BasepairStateOP const & {
		 return ptr.state(); } )
		.def("diff",[] (Basepair  & ptr, Basepair & bp) ->  float {
		 return ptr.diff(bp); } )
		.def("residues",[] (Basepair const & ptr) ->  ResidueOPs const {
		 return ptr.residues(); } )
		.def("partner",[] (Basepair  & ptr, ResidueOP const & res) ->  ResidueOP const & {
		 return ptr.partner(res); } )
		.def("name",[] (Basepair const & ptr) ->  String const {
		 return ptr.name(); } )
		.def("flip",[] (Basepair  & ptr) {
		ptr.flip(); } )
		.def("r",[] (Basepair const & ptr) ->  math::Matrix const & {
		 return ptr.r(); } )
		.def("d",[] (Basepair const & ptr) ->  math::Point const {
		 return ptr.d(); } )
		.def("uuid",[] (Basepair const & ptr) ->  util::Uuid const & {
		 return ptr.uuid(); } )
		.def("res1",[] (Basepair const & ptr) ->  ResidueOP {
		 return ptr.res1(); } )
		.def("res2",[] (Basepair const & ptr) ->  ResidueOP {
		 return ptr.res2(); } )
		.def("bp_type",[] (Basepair const & ptr) ->  String const & {
		 return ptr.bp_type(); } )
		.def("flipped",[] (Basepair const & ptr) ->  int const {
		 return ptr.flipped(); } )
		.def("atoms",[] (Basepair const & ptr) ->  AtomOPs const & {
		 return ptr.atoms(); } )
		.def("r",[] (Basepair  & ptr, math::Matrix const & nr) {
		ptr.r(nr); } )
		.def("uuid",[] (Basepair  & ptr, util::Uuid const & nuuid) {
		ptr.uuid(nuuid); } )
		.def("flipped",[] (Basepair  & ptr, int const & nflipped) {
		ptr.flipped(nflipped); } )
		.def("res1",[] (Basepair  & ptr, ResidueOP const & nres1) {
		ptr.res1(nres1); } )
		.def("res2",[] (Basepair  & ptr, ResidueOP const & nres2) {
		ptr.res2(nres2); } )
		.def("bp_type",[] (Basepair  & ptr, String const & nbp_type) {
		ptr.bp_type(nbp_type); } )
		.def("to_str",[] (Basepair const & ptr) -> String const {
		 return ptr.to_str(); } )
		.def("to_pdb_str",[] (Basepair const & ptr) -> String const {
		 return ptr.to_pdb_str(); } )
		.def("to_pdb",[] (Basepair const & ptr, String const fname) {
		ptr.to_pdb(fname); } )
		// operators
		.def(py::self == py::self)
		;

        py::class_<BasepairState, std::shared_ptr<BasepairState>>(m, "BasepairState")
		// ctors
		.def(py::init<>())
		.def(py::init<math::Point const &,math::Matrix const &,math::Points const &>())
		.def(py::init<BasepairState const &>())
		.def(py::init<String const &>())
		// methods
		.def("copy",[] (BasepairState const & ptr) ->  BasepairState {
		 return ptr.copy(); } )
		.def("move",[] (BasepairState  & ptr, math::Point const & p) {
		ptr.move(p); } )
		.def("transform",[] (BasepairState  & ptr, math::Matrix const & r, math::Vector const & t, math::Point & dummy) {
		ptr.transform(r, t, dummy); } )
		.def("transform",[] (BasepairState  & ptr, math::Matrix const & r, math::Vector const & t) {
		ptr.transform(r, t); } )
		.def("calculate_r_T",[] (BasepairState  & ptr) {
		ptr.calculate_r_T(); } )
		.def("get_transforming_r_and_t",[] (BasepairState  & ptr, BasepairState const & o_state, BasepairState & r_state) {
		ptr.get_transforming_r_and_t(o_state, r_state); } )
		.def("get_transformed_state",[] (BasepairState  & ptr, BasepairState const & o_state, BasepairState & r_state) {
		ptr.get_transformed_state(o_state, r_state); } )
		.def("flip",[] (BasepairState  & ptr) {
		ptr.flip(); } )
		.def("diff",[] (BasepairState  & ptr, BasepairStateOP const & state) ->  float {
		 return ptr.diff(state); } )
		.def("diff",[] (BasepairState  & ptr, BasepairState const & state) ->  float {
		 return ptr.diff(state); } )
		.def("_rot_diff",[] (BasepairState  & ptr, BasepairState const & state) ->  float {
		 return ptr._rot_diff(state); } )
		.def("to_str",[] (BasepairState const & ptr) ->  String const {
		 return ptr.to_str(); } )
		.def("d",[] (BasepairState  & ptr) ->  const math::Point & {
		 return ptr.d(); } )
		.def("d",[] (BasepairState const & ptr) ->  const math::Point & {
		 return ptr.d(); } )
		.def("r",[] (BasepairState  & ptr) ->  const math::Matrix & {
		 return ptr.r(); } )
		.def("r",[] (BasepairState const & ptr) ->  const math::Matrix & {
		 return ptr.r(); } )
		.def("r_T",[] (BasepairState  & ptr) ->  const math::Matrix & {
		 return ptr.r_T(); } )
		.def("r_T",[] (BasepairState const & ptr) ->  const math::Matrix & {
		 return ptr.r_T(); } )
		.def("sugars",[] (BasepairState  & ptr) ->  const math::Points & {
		 return ptr.sugars(); } )
		.def("sugars",[] (BasepairState const & ptr) ->  const math::Points & {
		 return ptr.sugars(); } )
		.def("d",[] (BasepairState  & ptr, math::Point const & newd) {
		ptr.d(newd); } )
		.def("r",[] (BasepairState  & ptr, math::Matrix const & newr) {
		ptr.r(newr); } )
		.def("sugars",[] (BasepairState  & ptr, math::Points const & newsug) {
		ptr.sugars(newsug); } )
		.def("set",[] (BasepairState  & ptr, BasepairState const & nstate_) {
		ptr.set(nstate_); } )
		;

        py::class_<Bead, std::shared_ptr<Bead>>(m, "Bead")
		// ctors
		.def(py::init<>())
		.def(py::init<math::Point const &,BeadType const>())
		.def(py::init<String const &>())
		.def(py::init<Bead const &>())
		// methods
		.def("distance",[] (Bead const & ptr, Bead const & b) ->  double {
		 return ptr.distance(b); } )
		.def("center",[] (Bead const & ptr) ->  math::Point {
		 return ptr.center(); } )
		.def("btype",[] (Bead const & ptr) ->  BeadType {
		 return ptr.btype(); } )
		;

        py::class_<BeadType, std::shared_ptr<BeadType>>(m, "BeadType")
		;

        py::class_<CIFParser, std::shared_ptr<CIFParser>>(m, "CIFParser")
		// ctors
		.def(py::init<>())
		// methods
		.def("parse",[] (CIFParser  & ptr, String const & pdb_file, int protein, int rna, int others) -> ResidueOPs const & {
		 return ptr.parse(pdb_file, protein, rna, others); } )
		;

        py::class_<Chain, std::shared_ptr<Chain>>(m, "Chain")
		// ctors
		.def(py::init<>())
		.def(py::init<ResidueOPs const &>())
		.def(py::init<Chain const &>())
		.def(py::init<String const &,ResidueTypeSet const &>())
		// methods
		.def("length",[] (Chain const & ptr) ->  int {
		 return ptr.length(); } )
		.def("first",[] (Chain  & ptr) ->  ResidueOP const & {
		 return ptr.first(); } )
		.def("last",[] (Chain  & ptr) ->  ResidueOP const & {
		 return ptr.last(); } )
		.def("subchain",[] (Chain  & ptr, int start, int end) ->  ChainOP {
		 return ptr.subchain(start, end); } )
		.def("subchain",[] (Chain  & ptr, ResidueOP const & r1, ResidueOP const & r2) ->  ChainOP {
		 return ptr.subchain(r1, r2); } )
		.def("to_str",[] (Chain const & ptr) -> String {
		 return ptr.to_str(); } )
		.def("to_pdb_str",[] (Chain const & ptr, int & acount, int rnum, String const & chain_id) -> String {
		 return ptr.to_pdb_str(acount, rnum, chain_id); } )
		.def("to_pdb_str",[] (Chain const & ptr, int & acount) ->  String {
		 return ptr.to_pdb_str(acount); } )
		.def("to_pdb",[] (Chain const & ptr, String const fname, int rnum, String const & chain_id) {
		ptr.to_pdb(fname, rnum, chain_id); } )
		.def("to_pdb",[] (Chain  & ptr, String const & fname) {
		ptr.to_pdb(fname); } )
		.def("residues",[] (Chain  & ptr) ->  ResidueOPs & {
		 return ptr.residues(); } )
		;

        py::class_<PDBParser, std::shared_ptr<PDBParser>>(m, "PDBParser")
		// ctors
		.def(py::init<>())
		// methods
		.def("parse",[] (PDBParser  & ptr, String const & pdb_file, int protein, int rna, int others) -> ResidueOPs const & {
		 return ptr.parse(pdb_file, protein, rna, others); } )
		;

        py::class_<RNAStructure, std::shared_ptr<RNAStructure>>(m, "RNAStructure")
		// ctors
		.def(py::init<>())
		.def(py::init<StructureOP const &,BasepairOPs const &,BasepairOPs const &>())
		.def(py::init<RNAStructure const &>())
		// methods
		.def("get_basepair",[] (RNAStructure  & ptr, String const & str) -> BasepairOPs {
		 return ptr.get_basepair(str); } )
		.def("get_basepair",[] (RNAStructure  & ptr, util::Uuid const & id) -> BasepairOPs {
		 return ptr.get_basepair(id); } )
		.def("get_basepair",[] (RNAStructure  & ptr, ResidueOP const & r1, ResidueOP const & r2) -> BasepairOPs {
		 return ptr.get_basepair(r1, r2); } )
		.def("get_basepair",[] (RNAStructure  & ptr, util::Uuid const & id1, util::Uuid const & id2) -> BasepairOPs {
		 return ptr.get_basepair(id1, id2); } )
		.def("get_beads",[] (RNAStructure  & ptr, BasepairOPs const & bps) -> Beads const & {
		 return ptr.get_beads(bps); } )
		.def("get_beads",[] (RNAStructure  & ptr, BasepairOP const & bps) -> Beads const & {
		 return ptr.get_beads(bps); } )
		.def("get_beads",[] (RNAStructure  & ptr) ->  Beads const & {
		 return ptr.get_beads(); } )
		.def("get_end_index",[] (RNAStructure  & ptr, BasepairOP const & bp) -> int {
		 return ptr.get_end_index(bp); } )
		.def("get_end_index",[] (RNAStructure  & ptr, String const & str) -> int {
		 return ptr.get_end_index(str); } )
		.def("to_pdb_str",[] (RNAStructure  & ptr, int rnumber, int close_chains) -> String const {
		 return ptr.to_pdb_str(rnumber, close_chains); } )
		.def("to_pdb",[] (RNAStructure  & ptr, String const fname, int renumber, int close_chains, int conect_statements) {
		ptr.to_pdb(fname, renumber, close_chains, conect_statements); } )
		.def("atoms",[] (RNAStructure const & ptr) ->  AtomOPs const {
		 return ptr.atoms(); } )
		.def("residues",[] (RNAStructure const & ptr) ->  ResidueOPs const {
		 return ptr.residues(); } )
		.def("chains",[] (RNAStructure  & ptr) -> ChainOPs const & {
		 return ptr.chains(); } )
		.def("get_residue",[] (RNAStructure  & ptr, int num, String const & chain_id, String const & i_code) ->  ResidueOP const {
		 return ptr.get_residue(num, chain_id, i_code); } )
		.def("get_residue",[] (RNAStructure  & ptr, util::Uuid const & uuid) ->  ResidueOP const {
		 return ptr.get_residue(uuid); } )
		.def("ends",[] (RNAStructure const & ptr) ->  BasepairOPs const & {
		 return ptr.ends(); } )
		.def("end_ids",[] (RNAStructure const & ptr) ->  Strings const & {
		 return ptr.end_ids(); } )
		.def("name",[] (RNAStructure const & ptr) ->  String const & {
		 return ptr.name(); } )
		.def("path",[] (RNAStructure const & ptr) ->  String const & {
		 return ptr.path(); } )
		.def("basepairs",[] (RNAStructure const & ptr) ->  BasepairOPs const & {
		 return ptr.basepairs(); } )
		.def("UNSAFE_basepairs",[] (RNAStructure  & ptr) ->  BasepairOPs & {
		 return ptr.UNSAFE_basepairs(); } )
		.def("beads",[] (RNAStructure const & ptr) ->  Beads const & {
		 return ptr.beads(); } )
		.def("score",[] (RNAStructure const & ptr) ->  float const & {
		 return ptr.score(); } )
		.def("protein_beads",[] (RNAStructure  & ptr) ->  Beads const & {
		 return ptr.protein_beads(); } )
		.def("name",[] (RNAStructure  & ptr, String const & nname) {
		ptr.name(nname); } )
		.def("path",[] (RNAStructure  & ptr, String const & npath) {
		ptr.path(npath); } )
		.def("score",[] (RNAStructure  & ptr, float const & nscore) {
		ptr.score(nscore); } )
		.def("end_ids",[] (RNAStructure  & ptr, Strings const & end_ids) {
		ptr.end_ids(end_ids); } )
		.def("ends",[] (RNAStructure  & ptr, BasepairOPs const & ends) {
		ptr.ends(ends); } )
		.def("protein_beads",[] (RNAStructure  & ptr, Beads const & beads) {
		ptr.protein_beads(beads); } )
		;

        py::class_<Residue, std::shared_ptr<Residue>>(m, "Residue")
		// ctors
		.def(py::init<ResidueType const &,String const &,int const &,String const &,String const &>())
		.def(py::init<Residue const &>())
		.def(py::init<String const &,ResidueTypeSet const &>())
		// methods
		.def("setup_atoms",[] (Residue  & ptr, AtomOPs & atoms) {
		ptr.setup_atoms(atoms); } )
		.def("get_atom",[] (Residue const & ptr, String const & name) ->  AtomOP const & {
		 return ptr.get_atom(name); } )
		.def("connected_to",[] (Residue const & ptr, Residue const & res, float cutoff) ->  int {
		 return ptr.connected_to(res, cutoff); } )
		.def("new_uuid",[] (Residue  & ptr) {
		ptr.new_uuid(); } )
		.def("get_beads",[] (Residue const & ptr) -> Beads {
		 return ptr.get_beads(); } )
		.def("to_str",[] (Residue const & ptr) -> String {
		 return ptr.to_str(); } )
		.def("to_pdb_str",[] (Residue  & ptr, int & acount) ->  String {
		 return ptr.to_pdb_str(acount); } )
		.def("to_pdb_str",[] (Residue const & ptr, int & acount, int rnum, String const & chain_id) -> String {
		 return ptr.to_pdb_str(acount, rnum, chain_id); } )
		.def("to_pdb",[] (Residue  & ptr, String const fname) {
		ptr.to_pdb(fname); } )
		.def("move",[] (Residue  & ptr, math::Point const & p) {
		ptr.move(p); } )
		.def("transform",[] (Residue  & ptr, math::Transform const & t) {
		ptr.transform(t); } )
		.def("num",[] (Residue  & ptr, int nnum) {
		ptr.num(nnum); } )
		.def("chain_id",[] (Residue  & ptr, String const & nchain_id) {
		ptr.chain_id(nchain_id); } )
		.def("uuid",[] (Residue  & ptr, util::Uuid const & uuid) {
		ptr.uuid(uuid); } )
		.def("center",[] (Residue const & ptr) ->  math::Point {
		 return ptr.center(); } )
		.def("name",[] (Residue const & ptr) ->  String const & {
		 return ptr.name(); } )
		.def("chain_id",[] (Residue const & ptr) ->  String const & {
		 return ptr.chain_id(); } )
		.def("i_code",[] (Residue const & ptr) ->  String const & {
		 return ptr.i_code(); } )
		.def("num",[] (Residue const & ptr) ->  int {
		 return ptr.num(); } )
		.def("short_name",[] (Residue const & ptr) ->  String {
		 return ptr.short_name(); } )
		.def("atoms",[] (Residue const & ptr) ->  AtomOPs const & {
		 return ptr.atoms(); } )
		.def("uuid",[] (Residue const & ptr) ->  util::Uuid const & {
		 return ptr.uuid(); } )
		// operators
		.def(py::self == py::self)
		;

        py::class_<ResidueType, std::shared_ptr<ResidueType>>(m, "ResidueType")
		// ctors
		.def(py::init<>())
		.def(py::init<String const,StringIntMap const,SetType const &>())
		// methods
		.def("get_correct_atom_name",[] (ResidueType const & ptr, Atom const & a) -> String {
		 return ptr.get_correct_atom_name(a); } )
		.def("match_name",[] (ResidueType const & ptr, String const & name) -> int {
		 return ptr.match_name(name); } )
		.def("short_name",[] (ResidueType const & ptr) ->  String {
		 return ptr.short_name(); } )
		.def("atom_pos_by_name",[] (ResidueType const & ptr, String const & aname) ->  int {
		 return ptr.atom_pos_by_name(aname); } )
		.def("size",[] (ResidueType  & ptr) ->  int {
		 return ptr.size(); } )
		.def("name",[] (ResidueType const & ptr) ->  String const & {
		 return ptr.name(); } )
		.def("set_type",[] (ResidueType const & ptr) ->  SetType {
		 return ptr.set_type(); } )
		;

        py::class_<ResidueTypeSet, std::shared_ptr<ResidueTypeSet>>(m, "ResidueTypeSet")
		// ctors
		.def(py::init<>())
		// methods
		.def("get_rtype_by_resname",[] (ResidueTypeSet const & ptr, String const & resname) -> ResidueType const & {
		 return ptr.get_rtype_by_resname(resname); } )
		.def("contains_rtype",[] (ResidueTypeSet const & ptr, String const & resname) -> bool  {
		 return ptr.contains_rtype(resname); } )
		;

        py::class_<ResidueTypeSetManager, std::unique_ptr<ResidueTypeSetManager,py::nodelete>>(m, "ResidueTypeSetManager")
	    .def(py::init([](){
	        return std::unique_ptr<ResidueTypeSetManager, py::nodelete>(&ResidueTypeSetManager::getInstance());
	    }))
        // methods
		.def("getInstance",[] (ResidueTypeSetManager  & ptr) -> ResidueTypeSetManager & {
		 return ptr.getInstance(); } )
		.def("residue_type_set",[] (ResidueTypeSetManager  & ptr) ->  ResidueTypeSet const & {
		 return ptr.residue_type_set(); } )
		;

        py::class_<SetType, std::shared_ptr<SetType>>(m, "SetType")
		;

        py::class_<Structure, std::shared_ptr<Structure>>(m, "Structure")
		// ctors
		.def(py::init<>())
		.def(py::init<ChainOPs const &>())
		.def(py::init<String const &>())
		.def(py::init<Structure const &>())
		.def(py::init<String const &,ResidueTypeSet const &>())
		// methods
		.def("renumber",[] (Structure  & ptr) {
		ptr.renumber(); } )
		.def("get_beads",[] (Structure  & ptr, ResidueOPs const & excluded_res) ->  Beads {
		 return ptr.get_beads(excluded_res); } )
		.def("get_beads",[] (Structure  & ptr) ->  Beads {
		 return ptr.get_beads(); } )
		.def("get_residue",[] (Structure  & ptr, int const & num, String const & chain_id, String const & i_code) -> ResidueOP {
		 return ptr.get_residue(num, chain_id, i_code); } )
		.def("get_residue",[] (Structure  & ptr, int const & num, String const & chain_id, String const & i_code) -> ResidueOP {
		 return ptr.get_residue(num, chain_id, i_code); } )
		.def("residues",[] (Structure const & ptr) -> ResidueOPs const {
		 return ptr.residues(); } )
		.def("atoms",[] (Structure  & ptr) ->  AtomOPs const {
		 return ptr.atoms(); } )
		.def("move",[] (Structure  & ptr, math::Point const & p) {
		ptr.move(p); } )
		.def("transform",[] (Structure  & ptr, math::Transform const & t) {
		ptr.transform(t); } )
		.def("to_pdb_str",[] (Structure  & ptr, int renumber, int conect_statements) -> String {
		 return ptr.to_pdb_str(renumber, conect_statements); } )
		.def("to_str",[] (Structure  & ptr) -> String {
		 return ptr.to_str(); } )
		.def("to_pdb",[] (Structure  & ptr, String const fname, int renumber, int conect_statements) {
		ptr.to_pdb(fname, renumber, conect_statements); } )
		.def("chains",[] (Structure const & ptr) ->  ChainOPs const & {
		 return ptr.chains(); } )
		;



}