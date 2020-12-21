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
using namespace io::error;

// Bindings start

namespace py = pybind11;

PYBIND11_MODULE(_io_error,m) {

	// classes
/*
        py::class_<base, std::shared_ptr<base>>(m, "base")
		// methods
		.def("format_error_message",[] (base const & ptr) {
		ptr.format_error_message(); } )
		.def("what",[] (base const & ptr) -> const char * {
		 return ptr.what(); } )
		// public attributes
		.def_readwrite("error_message_buffer", &base::error_message_buffer)
		;

        py::class_<can_not_open_file, std::shared_ptr<can_not_open_file>>(m, "can_not_open_file")
		// methods
		.def("format_error_message",[] (can_not_open_file const & ptr) {
		ptr.format_error_message(); } )
		// inherited methods
		.def("format_error_message",[] (base const & ptr) {
		ptr.format_error_message(); } )
		.def("what",[] (base const & ptr) -> const char * {
		 return ptr.what(); } )
		;

        py::class_<duplicated_column_in_header, std::shared_ptr<duplicated_column_in_header>>(m, "duplicated_column_in_header")
		// methods
		.def("format_error_message",[] (duplicated_column_in_header const & ptr) {
		ptr.format_error_message(); } )
		// inherited methods
		.def("format_error_message",[] (base const & ptr) {
		ptr.format_error_message(); } )
		.def("what",[] (base const & ptr) -> const char * {
		 return ptr.what(); } )
		;

        py::class_<escaped_string_not_closed, std::shared_ptr<escaped_string_not_closed>>(m, "escaped_string_not_closed")
		// methods
		.def("format_error_message",[] (escaped_string_not_closed const & ptr) {
		ptr.format_error_message(); } )
		// inherited methods
		.def("format_error_message",[] (base const & ptr) {
		ptr.format_error_message(); } )
		.def("what",[] (base const & ptr) -> const char * {
		 return ptr.what(); } )
		;

        py::class_<extra_column_in_header, std::shared_ptr<extra_column_in_header>>(m, "extra_column_in_header")
		// methods
		.def("format_error_message",[] (extra_column_in_header const & ptr) {
		ptr.format_error_message(); } )
		// inherited methods
		.def("format_error_message",[] (base const & ptr) {
		ptr.format_error_message(); } )
		.def("what",[] (base const & ptr) -> const char * {
		 return ptr.what(); } )
		;

        py::class_<header_missing, std::shared_ptr<header_missing>>(m, "header_missing")
		// methods
		.def("format_error_message",[] (header_missing const & ptr) {
		ptr.format_error_message(); } )
		// inherited methods
		.def("format_error_message",[] (base const & ptr) {
		ptr.format_error_message(); } )
		.def("what",[] (base const & ptr) -> const char * {
		 return ptr.what(); } )
		;

        py::class_<integer_must_be_positive, std::shared_ptr<integer_must_be_positive>>(m, "integer_must_be_positive")
		// methods
		.def("format_error_message",[] (integer_must_be_positive const & ptr) {
		ptr.format_error_message(); } )
		// inherited methods
		.def("format_error_message",[] (base const & ptr) {
		ptr.format_error_message(); } )
		.def("what",[] (base const & ptr) -> const char * {
		 return ptr.what(); } )
		;

        py::class_<integer_overflow, std::shared_ptr<integer_overflow>>(m, "integer_overflow")
		// methods
		.def("format_error_message",[] (integer_overflow const & ptr) {
		ptr.format_error_message(); } )
		// inherited methods
		.def("format_error_message",[] (base const & ptr) {
		ptr.format_error_message(); } )
		.def("what",[] (base const & ptr) -> const char * {
		 return ptr.what(); } )
		;

        py::class_<integer_underflow, std::shared_ptr<integer_underflow>>(m, "integer_underflow")
		// methods
		.def("format_error_message",[] (integer_underflow const & ptr) {
		ptr.format_error_message(); } )
		// inherited methods
		.def("format_error_message",[] (base const & ptr) {
		ptr.format_error_message(); } )
		.def("what",[] (base const & ptr) -> const char * {
		 return ptr.what(); } )
		;

        py::class_<invalid_single_character, std::shared_ptr<invalid_single_character>>(m, "invalid_single_character")
		// methods
		.def("format_error_message",[] (invalid_single_character const & ptr) {
		ptr.format_error_message(); } )
		// inherited methods
		.def("format_error_message",[] (base const & ptr) {
		ptr.format_error_message(); } )
		.def("what",[] (base const & ptr) -> const char * {
		 return ptr.what(); } )
		;

        py::class_<line_length_limit_exceeded, std::shared_ptr<line_length_limit_exceeded>>(m, "line_length_limit_exceeded")
		// methods
		.def("format_error_message",[] (line_length_limit_exceeded const & ptr) {
		ptr.format_error_message(); } )
		// inherited methods
		.def("format_error_message",[] (base const & ptr) {
		ptr.format_error_message(); } )
		.def("what",[] (base const & ptr) -> const char * {
		 return ptr.what(); } )
		;

        py::class_<missing_column_in_header, std::shared_ptr<missing_column_in_header>>(m, "missing_column_in_header")
		// methods
		.def("format_error_message",[] (missing_column_in_header const & ptr) {
		ptr.format_error_message(); } )
		// inherited methods
		.def("format_error_message",[] (base const & ptr) {
		ptr.format_error_message(); } )
		.def("what",[] (base const & ptr) -> const char * {
		 return ptr.what(); } )
		;

        py::class_<no_digit, std::shared_ptr<no_digit>>(m, "no_digit")
		// methods
		.def("format_error_message",[] (no_digit const & ptr) {
		ptr.format_error_message(); } )
		// inherited methods
		.def("format_error_message",[] (base const & ptr) {
		ptr.format_error_message(); } )
		.def("what",[] (base const & ptr) -> const char * {
		 return ptr.what(); } )
		;

        py::class_<too_few_columns, std::shared_ptr<too_few_columns>>(m, "too_few_columns")
		// methods
		.def("format_error_message",[] (too_few_columns const & ptr) {
		ptr.format_error_message(); } )
		// inherited methods
		.def("format_error_message",[] (base const & ptr) {
		ptr.format_error_message(); } )
		.def("what",[] (base const & ptr) -> const char * {
		 return ptr.what(); } )
		;

        py::class_<too_many_columns, std::shared_ptr<too_many_columns>>(m, "too_many_columns")
		// methods
		.def("format_error_message",[] (too_many_columns const & ptr) {
		ptr.format_error_message(); } )
		// inherited methods
		.def("format_error_message",[] (base const & ptr) {
		ptr.format_error_message(); } )
		.def("what",[] (base const & ptr) -> const char * {
		 return ptr.what(); } )
		;

        py::class_<with_column_content, std::shared_ptr<with_column_content>>(m, "with_column_content")
		// ctors
		.def(py::init<>())
		// methods
		.def("set_column_content",[] (with_column_content  & ptr, const char * column_content) {
		ptr.set_column_content(column_content); } )
		// public attributes
		.def_readwrite("column_content", &with_column_content::column_content)
		;

        py::class_<with_column_name, std::shared_ptr<with_column_name>>(m, "with_column_name")
		// ctors
		.def(py::init<>())
		// methods
		.def("set_column_name",[] (with_column_name  & ptr, const char * column_name) {
		ptr.set_column_name(column_name); } )
		// public attributes
		.def_readwrite("column_name", &with_column_name::column_name)
		;

        py::class_<with_errno, std::shared_ptr<with_errno>>(m, "with_errno")
		// ctors
		.def(py::init<>())
		// methods
		.def("set_errno",[] (with_errno  & ptr, int errno_value) {
		ptr.set_errno(errno_value); } )
		// public attributes
		.def_readwrite("errno_value", &with_errno::errno_value)
		;

        py::class_<with_file_line, std::shared_ptr<with_file_line>>(m, "with_file_line")
		// ctors
		.def(py::init<>())
		// methods
		.def("set_file_line",[] (with_file_line  & ptr, int file_line) {
		ptr.set_file_line(file_line); } )
		// public attributes
		.def_readwrite("file_line", &with_file_line::file_line)
		;

        py::class_<with_file_name, std::shared_ptr<with_file_name>>(m, "with_file_name")
		// ctors
		.def(py::init<>())
		// methods
		.def("set_file_name",[] (with_file_name  & ptr, const char * file_name) {
		ptr.set_file_name(file_name); } )
		// public attributes
		.def_readwrite("file_name", &with_file_name::file_name)
		;

		*/


}