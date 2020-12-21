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
using namespace util;

// Bindings start

namespace py = pybind11;

PYBIND11_MODULE(_util,m) {

	// free functions

    m.def("compare_bps", [] (X3dna::X3Basepairs & lhs, X3dna::X3Basepairs & rhs) -> String {
		 return compare_bps(lhs, rhs); },
	py::arg("lhs"),
	py::arg("rhs") 
    );

    m.def("get_str_from_x3dna_type", [] (X3dnaBPType type) -> String {
		 return get_str_from_x3dna_type(type); },
	py::arg("type") 
    );

    m.def("get_x3dna_by_type", [] (String const & name) -> X3dnaBPType {
		 return get_x3dna_by_type(name); },
	py::arg("name") 
    );

    m.def("json_cleanup", [] () {
		json_cleanup(); } 
    );

    m.def("points_to_pdb", [] (String const & filename, math::Points const & points) {
		points_to_pdb(filename, points); },
	py::arg("filename"),
	py::arg("points") 
    );

    m.def("points_to_pdb_str", [] (math::Points const & points) -> String {
		 return points_to_pdb_str(points); },
	py::arg("points") 
    );

    m.def("str_to_type", [] (String const s) -> MotifType const {
		 return str_to_type(s); },
	py::arg("s") 
    );

    m.def("type_to_str", [] (MotifType const mtype) -> String const {
		 return type_to_str(mtype); },
	py::arg("mtype") 
    );
	// classes
#define Value double
#define Values std::vector<Value>
        py::class_<CartesianProduct<Value>, std::shared_ptr<CartesianProduct<Value>>>(m, "CartesianProduct<Value>")
		// ctors
		.def(py::init<>())
		.def(py::init<std::vector<Values> const &>())
		// methods
		.def("setup",[] (CartesianProduct<Value>  & ptr, std::vector<Values> const & values) {
		ptr.setup(values); } )
		.def("end",[] (CartesianProduct<Value>  & ptr) ->  int const {
		 return ptr.end(); } )
		.def("next",[] (CartesianProduct<Value>  & ptr) -> Values const & {
		 return ptr.next(); } )
		;
#undef Values
#undef Value
        py::class_<MonteCarlo, std::shared_ptr<MonteCarlo>>(m, "MonteCarlo")
		// ctors
		.def(py::init<float>())
		// methods
		.def("accept",[] (MonteCarlo  & ptr, float current, float next) ->  int {
		 return ptr.accept(current, next); } )
		.def("set_temperature",[] (MonteCarlo  & ptr, float new_temp) {
		ptr.set_temperature(new_temp); } )
		.def("get_temperature",[] (MonteCarlo  & ptr) ->  float {
		 return ptr.get_temperature(); } )
		.def("scale_temperature",[] (MonteCarlo  & ptr, float scale) {
		ptr.scale_temperature(scale); } )
		;

        py::class_<RandomNumberGenerator, std::shared_ptr<RandomNumberGenerator>>(m, "RandomNumberGenerator")
		// ctors
		.def(py::init<>())
		// methods
		.def("rand",[] (RandomNumberGenerator  & ptr) ->  float {
		 return ptr.rand(); } )
		.def("randrange",[] (RandomNumberGenerator  & ptr, int i) ->  int {
		 return ptr.randrange(i); } )
		;

        py::class_<Sqlite3Connection, std::shared_ptr<Sqlite3Connection>>(m, "Sqlite3Connection")
		// ctors
		.def(py::init<>())
		.def(py::init<String const>())
		// methods
		.def("query",[] (Sqlite3Connection  & ptr, String const & query_statement) {
		ptr.query(query_statement); } )
		.def("count",[] (Sqlite3Connection  & ptr) -> int {
		 return ptr.count(); } )
		.def("fetch_one",[] (Sqlite3Connection  & ptr, String const & query_statement) -> Strings {
		 return ptr.fetch_one(query_statement); } )
		;

        py::class_<StericLookup, std::shared_ptr<StericLookup>>(m, "StericLookup")
		// ctors
		.def(py::init<>())
		.def(py::init<float,float,int>())
		// methods
		.def("add_point",[] (StericLookup  & ptr, math::Point const & p) {
		ptr.add_point(p); } )
		.def("add_points",[] (StericLookup  & ptr, math::Points const & points) {
		ptr.add_points(points); } )
		.def("clash",[] (StericLookup  & ptr, math::Point const & p) -> int {
		 return ptr.clash(p); } )
		.def("clash",[] (StericLookup  & ptr, math::Point const & p) -> int {
		 return ptr.clash(p); } )
		.def("better_clash",[] (StericLookup  & ptr, math::Point const & p) -> int {
		 return ptr.better_clash(p); } )
		.def("total_clash",[] (StericLookup  & ptr, math::Point const & p) -> int {
		 return ptr.total_clash(p); } )
		.def("total_clash",[] (StericLookup  & ptr, math::Point const & p) -> int {
		 return ptr.total_clash(p); } )
		;

        py::class_<StericLookupNew, std::shared_ptr<StericLookupNew>>(m, "StericLookupNew")
		// ctors
		.def(py::init<>())
		// methods
		.def("add_point",[] (StericLookupNew  & ptr, math::Point const & p) {
		ptr.add_point(p); } )
		.def("add_points",[] (StericLookupNew  & ptr, math::Points const & points) {
		ptr.add_points(points); } )
		.def("clash",[] (StericLookupNew  & ptr, math::Point const & p) -> bool {
		 return ptr.clash(p); } )
		.def("clash",[] (StericLookupNew  & ptr, math::Point const & p) -> bool {
		 return ptr.clash(p); } )
		.def("to_pdb",[] (StericLookupNew  & ptr, String const & pdb_name) {
		ptr.to_pdb(pdb_name); } )
		.def("size",[] (StericLookupNew  & ptr) -> int {
		 return ptr.size(); } )
		;

        py::class_<Uuid, std::shared_ptr<Uuid>>(m, "Uuid")
		// ctors
		.def(py::init<>())
		// methods
		.def("s_uuid",[] (Uuid const & ptr) ->  String const & {
		 return ptr.s_uuid(); } )
		// operators
		.def(py::self == py::self)
		.def(py::self == py::self)
		.def(py::self != py::self)
		.def(py::self != py::self)
		.def(py::self < py::self)
		;
/*
        py::class_<UuidCompare, std::shared_ptr<UuidCompare>>(m, "UuidCompare")
		// operators
		.def(py::self () Uuid const &)
		;
*/
        py::class_<X3dna::X3BPInfo, std::shared_ptr<X3dna::X3BPInfo>>(m, "X3BPInfo")
		// public attributes
		.def_readwrite("res1_type_name", &X3dna::X3BPInfo::res1_type_name)
		.def_readwrite("res2_type_name", &X3dna::X3BPInfo::res2_type_name)
		.def_readwrite("res1_name", &X3dna::X3BPInfo::res1_name)
		.def_readwrite("res2_name", &X3dna::X3BPInfo::res2_name)
		.def_readwrite("res1_chain_id", &X3dna::X3BPInfo::res1_chain_id)
		.def_readwrite("res2_chain_id", &X3dna::X3BPInfo::res2_chain_id)
		.def_readwrite("res1_num", &X3dna::X3BPInfo::res1_num)
		.def_readwrite("res2_num", &X3dna::X3BPInfo::res2_num)
		;

        py::class_<X3dna::X3Basepair, std::shared_ptr<X3dna::X3Basepair>>(m, "X3Basepair")
		// methods
		.def("valid",[] (X3dna::X3Basepair const & ptr) -> bool {
		 return ptr.valid(); } )
		.def("to_string",[] (X3dna::X3Basepair const & ptr) -> String {
		 return ptr.to_string(); } )
		.def("key",[] (X3dna::X3Basepair const & ptr) -> unsigned int {
		 return ptr.key(); } )
		// public attributes
		.def_readwrite("res1", &X3dna::X3Basepair::res1)
		.def_readwrite("res2", &X3dna::X3Basepair::res2)
		.def_readwrite("d", &X3dna::X3Basepair::d)
		.def_readwrite("r", &X3dna::X3Basepair::r)
		.def_readwrite("bp_type", &X3dna::X3Basepair::bp_type)
		;

        py::class_<X3dna::X3Motif, std::shared_ptr<X3dna::X3Motif>>(m, "X3Motif")
		// public attributes
		.def_readwrite("residues", &X3dna::X3Motif::residues)
		.def_readwrite("mtype", &X3dna::X3Motif::mtype)
		;

        py::class_<X3dna::X3Residue, std::shared_ptr<X3dna::X3Residue>>(m, "X3Residue")
		// methods
		.def("valid",[] (X3dna::X3Residue const & ptr) -> bool {
		 return ptr.valid(); } )
		// operators
		.def(py::self == py::self)
		// public attributes
		.def_readwrite("num", &X3dna::X3Residue::num)
		.def_readwrite("chain_id", &X3dna::X3Residue::chain_id)
		.def_readwrite("i_code", &X3dna::X3Residue::i_code)
		;

        py::class_<X3dna, std::shared_ptr<X3dna>>(m, "X3dna")
		// ctors
		.def(py::init<>())
		// methods
		.def("set_envs",[] (X3dna  & ptr) {
		ptr.set_envs(); } )
		.def("get_basepairs",[] (X3dna const & ptr, String const & pdb_path) -> X3dna::X3Basepairs {
		 return ptr.get_basepairs(pdb_path); } )
		.def("get_basepairs_json",[] (X3dna const & ptr, String const & pdb_path) -> X3dna::X3Basepairs  {
		 return ptr.get_basepairs_json(pdb_path); } )
		.def("get_motifs",[] (X3dna const & ptr, String const & strings) -> X3dna::X3Motifs {
		 return ptr.get_motifs(strings); } )
		.def("get_motifs",[] (X3dna const & ptr, String const & str) -> X3dna::X3Motifs {
		 return ptr.get_motifs(str); } )
		.def("set_rebuild_files",[] (X3dna const & ptr, bool rebuild_files) {
		ptr.set_rebuild_files(rebuild_files); } )
		;

	// enums
	py::enum_<MotifType>(m,"MotifType")
		.value("TWOWAY",MotifType::TWOWAY)
		.value("NWAY",MotifType::NWAY)
		.value("HAIRPIN",MotifType::HAIRPIN)
		.value("TCONTACT_HP_HP",MotifType::TCONTACT_HP_HP)
		.value("TCONTACT_H_HP",MotifType::TCONTACT_H_HP)
		.value("TCONTACT_H_H",MotifType::TCONTACT_H_H)
		.value("T_T",MotifType::T_T)
		.value("T_T_T",MotifType::T_T_T)
		.value("TWOWAY_SEGMENTS",MotifType::TWOWAY_SEGMENTS)
		.value("HELIX",MotifType::HELIX)
		.value("SSTRAND",MotifType::SSTRAND)
		.value("TCONTACT",MotifType::TCONTACT)
		.value("UNKNOWN",MotifType::UNKNOWN)
		.value("ALL",MotifType::ALL)
	;
	py::enum_<X3dnaBPType>(m,"X3dnaBPType")
		.value("cmU",X3dnaBPType::cmU)
		.value("cMUM",X3dnaBPType::cMUM)
		.value("tWPW",X3dnaBPType::tWPW)
		.value("cDPM",X3dnaBPType::cDPM)
		.value("DWPW",X3dnaBPType::DWPW)
		.value("tWUM",X3dnaBPType::tWUM)
		.value("tmUM",X3dnaBPType::tmUM)
		.value("cWPM",X3dnaBPType::cWPM)
		.value("DWUW",X3dnaBPType::DWUW)
		.value("cMPD",X3dnaBPType::cMPD)
		.value("cDUm",X3dnaBPType::cDUm)
		.value("cMPW",X3dnaBPType::cMPW)
		.value("tMPm",X3dnaBPType::tMPm)
		.value("tMUW",X3dnaBPType::tMUW)
		.value("cmUm",X3dnaBPType::cmUm)
		.value("cMUW",X3dnaBPType::cMUW)
		.value("cWUW",X3dnaBPType::cWUW)
		.value("cDUM",X3dnaBPType::cDUM)
		.value("cmPM",X3dnaBPType::cmPM)
		.value("cmUM",X3dnaBPType::cmUM)
		.value("DDDD",X3dnaBPType::DDDD)
		.value("cmUW",X3dnaBPType::cmUW)
		.value("tMUm",X3dnaBPType::tMUm)
		.value("cDUW",X3dnaBPType::cDUW)
		.value("cMPm",X3dnaBPType::cMPm)
		.value("cMUm",X3dnaBPType::cMUm)
		.value("cDDD",X3dnaBPType::cDDD)
		.value("tWPm",X3dnaBPType::tWPm)
		.value("cDPm",X3dnaBPType::cDPm)
		.value("tmPm",X3dnaBPType::tmPm)
		.value("tWPD",X3dnaBPType::tWPD)
		.value("tmPW",X3dnaBPType::tmPW)
		.value("tDDD",X3dnaBPType::tDDD)
		.value("cWUD",X3dnaBPType::cWUD)
		.value("cWUM",X3dnaBPType::cWUM)
		.value("tDUW",X3dnaBPType::tDUW)
		.value("tMPM",X3dnaBPType::tMPM)
		.value("tDUM",X3dnaBPType::tDUM)
		.value("cMUD",X3dnaBPType::cMUD)
		.value("cWUm",X3dnaBPType::cWUm)
		.value("tDPm",X3dnaBPType::tDPm)
		.value("tMUD",X3dnaBPType::tMUD)
		.value("cmPW",X3dnaBPType::cmPW)
		.value("cMPM",X3dnaBPType::cMPM)
		.value("cmPD",X3dnaBPType::cmPD)
		.value("cmUD",X3dnaBPType::cmUD)
		.value("cDUD",X3dnaBPType::cDUD)
		.value("cWPW",X3dnaBPType::cWPW)
		.value("tDUD",X3dnaBPType::tDUD)
		.value("tDPW",X3dnaBPType::tDPW)
		.value("tmUm",X3dnaBPType::tmUm)
		.value("cWPD",X3dnaBPType::cWPD)
		.value("tmPD",X3dnaBPType::tmPD)
		.value("tDPD",X3dnaBPType::tDPD)
		.value("cDPD",X3dnaBPType::cDPD)
		.value("tDUm",X3dnaBPType::tDUm)
		.value("tDPM",X3dnaBPType::tDPM)
		.value("tWUD",X3dnaBPType::tWUD)
		.value("tmUW",X3dnaBPType::tmUW)
		.value("tMUM",X3dnaBPType::tMUM)
		.value("tMPD",X3dnaBPType::tMPD)
		.value("cDPW",X3dnaBPType::cDPW)
		.value("tmPM",X3dnaBPType::tmPM)
		.value("tWUm",X3dnaBPType::tWUm)
		.value("cWPm",X3dnaBPType::cWPm)
		.value("tmUD",X3dnaBPType::tmUD)
		.value("tWPM",X3dnaBPType::tWPM)
		.value("DWPm",X3dnaBPType::DWPm)
		.value("tMPW",X3dnaBPType::tMPW)
		.value("DDPm",X3dnaBPType::DDPm)
		.value("tWUW",X3dnaBPType::tWUW)
		.value("cmPm",X3dnaBPType::cmPm)
		.value("DWUm",X3dnaBPType::DWUm)
		.value("DMPm",X3dnaBPType::DMPm)
		.value("DWPM",X3dnaBPType::DWPM)
		.value("DMPM",X3dnaBPType::DMPM)
		.value("DmPW",X3dnaBPType::DmPW)
		.value("DWUM",X3dnaBPType::DWUM)
		.value("DmPm",X3dnaBPType::DmPm)
		.value("DDUM",X3dnaBPType::DDUM)
		.value("DMUm",X3dnaBPType::DMUm)
		.value("DDUm",X3dnaBPType::DDUm)
		.value("DMPW",X3dnaBPType::DMPW)
		.value("DMPD",X3dnaBPType::DMPD)
		.value("DMUM",X3dnaBPType::DMUM)
		.value("DmUm",X3dnaBPType::DmUm)
		.value("DMUW",X3dnaBPType::DMUW)
		.value("DWUD",X3dnaBPType::DWUD)
	;


}