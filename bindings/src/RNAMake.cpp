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


// namespaces
using namespace base;
using namespace math;
using namespace data_structure;
using namespace util;
using namespace secondary_structure;
using namespace eternabot;
using namespace structure;
using namespace motif;
using namespace motif_tools;
using namespace resources;
using namespace motif_data_structure;
using namespace motif_search;
using namespace motif_search::path_finding;
using namespace thermo_fluctuation;

// Bindings start

namespace py = pybind11;

PYBIND11_MODULE(RNAMake,m) {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// base
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // free functions

    m.def("base_dir", [] (String const & path) -> String {
              return base_dir(path); },
          py::arg("path")
    );

    m.def("demangle", [] (std::string & str) -> std::string {
              return demangle(str); },
          py::arg("str")
    );

    m.def("determine_string_data_type", [] (String const & str) -> DataType {
              return determine_string_data_type(str); },
          py::arg("str")
    );

    m.def("execute_command", [] (const char * cmd) -> String {
              return execute_command(cmd); },
          py::arg("cmd")
    );

    m.def("execute_command_json", [] (const char * cmd) -> nlohmann::json {
              return execute_command_json(cmd); },
          py::arg("cmd")
    );

    m.def("file_exists", [] (String const & name) ->  bool {
              return file_exists(name); },
          py::arg("name")
    );

    m.def("filename", [] (String const & path) -> String {
              return filename(path); },
          py::arg("path")
    );

    m.def("get_lines_from_file", [] (String const& name)  -> Strings {
              return get_lines_from_file(name); },
          py::arg("noexcept(false")
    );

    m.def("get_os_name", [] () -> String {
        return get_os_name(); }
    );

    py::enum_<LogLevel>(m, "LogLevel")
            .value("DEBUG", LogLevel::DEBUG)
            .value("ERROR", LogLevel::ERROR)
            .value("FATAL", LogLevel::FATAL)
            .value("INFO", LogLevel::INFO)
            .value("VERBOSE", LogLevel::VERBOSE)
            .value("WARN", LogLevel::WARN)
            ;

    m.def("init_logging", [] (LogLevel log_level) { // TODO add the log level stuff
              init_logging(log_level); },
          py::arg("log_level") = LogLevel::INFO
    );

    m.def("is_dir", [] (String const & path) ->  int {
              return is_dir(path); },
          py::arg("path")
    );

    m.def("is_number", [] (String const & s) -> bool {
              return is_number(s); },
          py::arg("s")
    );

    m.def("join_by_delimiter", [] (Strings const & strs, String const & delimiter) -> String {
              return join_by_delimiter(strs, delimiter); },
          py::arg("strs"),
          py::arg("delimiter")
    );

    m.def("lib_path", [] () -> String {
        return lib_path(); }
    );

    m.def("log_level_from_str", [] (String const & s) -> LogLevel {
              return log_level_from_str(s); },
          py::arg("s")
    );

    m.def("ltrim", [] (String & s) -> String & {
              return ltrim(s); },
          py::arg("s")
    );

    m.def("motif_dirs", [] () -> String {
        return motif_dirs(); }
    );

    m.def("print_backtrace", [] () {
        print_backtrace(); }
    );

    m.def("replace_all", [] (String & context, String const & from, String const & to) -> String & {
              return replace_all(context, from, to); },
          py::arg("context"),
          py::arg("from"),
          py::arg("to")
    );

    m.def("resources_path", [] () -> String {
        return resources_path(); }
    );

    m.def("rtrim", [] (String & s) -> String & {
              return rtrim(s); },
          py::arg("s")
    );

    m.def("save_backtrace", [] () {
        save_backtrace(); }
    );

    m.def("split_str_by_delimiter", [] (String s, String delimiter) -> Strings {
              return split_str_by_delimiter(s, delimiter); },
          py::arg("s"),
          py::arg("delimiter")
    );

    m.def("tokenize_line", [] (String const & raw_line) -> Strings {
              return tokenize_line(raw_line); },
          py::arg("raw_line")
    );

    m.def("trim", [] (String & s) -> String & {
              return trim(s); },
          py::arg("s")
    );

    m.def("unittest_resource_dir", [] () -> String {
        return unittest_resource_dir(); }
    );

    m.def("x3dna_path", [] () -> String {
        return x3dna_path(); }
    );
    // classes
/*
        py::class_<Application, std::shared_ptr<Application>>(m, "Application")
		// ctors
		.def(py::init<>())
		// methods
		.def("setup_options",[] (Application  & ptr) {
		ptr.setup_options(); } )
		.def("parse_command_line",[] (Application  & ptr, int argc, const char * * argv) {
		ptr.parse_command_line(argc, argv); } )
		.def("run",[] (Application  & ptr) {
		ptr.run(); } )
		.def("add_option_int",[] (Application  & ptr, String const & name, int const & val, OptionType const & type, bool required) { ptr.add_option(name, val, type, required); } )
		.def("add_option_float",[] (Application  & ptr, String const & name, float const & val, OptionType const & type, bool required) { ptr.add_option(name, val, type, required); } )
		.def("add_option_bool",[] (Application  & ptr, String const & name, bool const & val, OptionType const & type, bool required) { ptr.add_option(name, val, type, required); } )
		.def("add_option_string",[] (Application  & ptr, String const & name, String const & val, OptionType const & type, bool required) { ptr.add_option(name, val, type, required); } )
		.def("add_cl_options",[] (Application  & ptr, Options const & opts, String prefix) {
		ptr.add_cl_options(opts, prefix); } )
		.def("get_int_option",[] (Application  & ptr, String const & name) ->  float {
		 return ptr.get_int_option(name); } )
		.def("get_float_option",[] (Application  & ptr, String const & name) ->  float {
		 return ptr.get_float_option(name); } )
		.def("get_string_option",[] (Application  & ptr, String const & name) ->  String {
		 return ptr.get_string_option(name); } )
		.def("get_bool_option",[] (Application  & ptr, String const & name) ->  bool {
		 return ptr.get_bool_option(name); } )
		;
*/
    py::class_<CommandLineOption, std::shared_ptr<CommandLineOption>>(m, "CommandLineOption")
            // ctors
            .def(py::init<String const &, int const &,OptionType const &,bool>())
            .def(py::init<String const &, float const &,OptionType const &,bool>())
            .def(py::init<String const &, bool const &,OptionType const &,bool>())
            .def(py::init<String const &, String const &,OptionType const &,bool>())
            .def(py::init<Option const &>())
                    // methods
            .def("filled",[] (CommandLineOption  & ptr) ->  bool {
                return ptr.filled(); } )
            .def("required",[] (CommandLineOption  & ptr) ->  bool {
                return ptr.required(); } )
            .def("filled",[] (CommandLineOption  & ptr, bool const & filled) {
                ptr.filled(filled); } )
                    // inherited methods
            .def("get_float",[] (Option  & ptr) -> float {
                return ptr.get_float(); } )
            .def("get_int",[] (Option  & ptr) -> int {
                return ptr.get_int(); } )
            .def("get_string",[] (Option  & ptr) -> String const & {
                return ptr.get_string(); } )
            .def("get_bool",[] (Option  & ptr) -> bool const & {
                return ptr.get_bool(); } )
            .def("value",[] (Option  & ptr, int const & i) {
                ptr.value(i); } )
            .def("value",[] (Option  & ptr, char const * c_str) {
                ptr.value(c_str); } )
            .def("value",[] (Option  & ptr, int const & i) {
                ptr.value(i); } )
            .def("value",[] (Option  & ptr, int const & i) {
                ptr.value(i); } )
            .def("value",[] (Option  & ptr, int const & i) {
                ptr.value(i); } )
            .def("name",[] (Option const & ptr) ->  String const & {
                return ptr.name(); } )
            .def("type",[] (Option const & ptr) ->  OptionType const & {
                return ptr.type(); } )
            .def("type_name",[] (Option  & ptr) ->  String {
                return ptr.type_name(); } )
            ;

    py::class_<CommandLineOptions, std::shared_ptr<CommandLineOptions>>(m, "CommandLineOptions")
            // ctors
            .def(py::init<>())
                    // methods
            .def("begin",[] (CommandLineOptions  & ptr) -> std::vector<CommandLineOptionOP>::iterator {
                return ptr.begin(); } )
            .def("end",[] (CommandLineOptions  & ptr) -> std::vector<CommandLineOptionOP>::iterator {
                return ptr.end(); } )
            .def("begin",[] (CommandLineOptions const & ptr) -> std::vector<CommandLineOptionOP>::const_iterator {
                return ptr.begin(); } )
            .def("end",[] (CommandLineOptions const & ptr) -> std::vector<CommandLineOptionOP>::const_iterator {
                return ptr.end(); } )
            .def("add_option",[] (CommandLineOptions  & ptr, String const & name, String const & value, OptionType const & type, bool required) { ptr.add_option(name, value, type, required); } )
            .def("add_option",[] (CommandLineOptions  & ptr, String const & name, int const & value, OptionType const & type, bool required) { ptr.add_option(name, value, type, required); } )
            .def("add_option",[] (CommandLineOptions  & ptr, String const & name, bool const & value, OptionType const & type, bool required) { ptr.add_option(name, value, type, required); } )
            .def("add_option",[] (CommandLineOptions  & ptr, String const & name, float const & value, OptionType const & type, bool required) { ptr.add_option(name, value, type, required); } )
            .def("add_options",[] (CommandLineOptions  & ptr, Options const & opts) {
                ptr.add_options(opts); } )
            .def("parse_command_line", [] (CommandLineOptions  & ptr, const int  argc, Strings & vec) {
                std::vector<char*> ptrs; ptrs.reserve(vec.size()); for(auto& v : vec) {ptrs.push_back(const_cast<char*>(v.c_str()));}
                ; ptr.parse_command_line(argc, (const char**)(&ptrs[0])); } )
            .def("get_int",[] (CommandLineOptions const & ptr, String const & name) ->  float {
                return ptr.get_int(name); } )
            .def("get_float",[] (CommandLineOptions const & ptr, String const & name) ->  float {
                return ptr.get_float(name); } )
            .def("get_string",[] (CommandLineOptions const & ptr, String const & name) ->  String {
                return ptr.get_string(name); } )
            .def("get_bool",[] (CommandLineOptions const & ptr, String const & name) ->  bool {
                return ptr.get_bool(name); } )
            .def("has_option",[] (CommandLineOptions const & ptr, String const & name) ->  bool {
                return ptr.has_option(name); } )
            .def("set_value",[] (CommandLineOptions  & ptr, String const & name, String const & val) { ptr.set_value(name, val); } )
            .def("set_value",[] (CommandLineOptions  & ptr, String const & name, int const & val) { ptr.set_value(name, val); } )
            .def("set_value",[] (CommandLineOptions  & ptr, String const & name, bool const & val) { ptr.set_value(name, val); } )
            .def("set_value",[] (CommandLineOptions  & ptr, String const & name, float const & val) { ptr.set_value(name, val); } )
            .def("is_filled",[] (CommandLineOptions const & ptr, String const & name) ->  bool {
                return ptr.is_filled(name); } )
            ;

    py::class_<CommandLineParser, std::shared_ptr<CommandLineParser>>(m, "CommandLineParser")
            // ctors
            .def(py::init<>())
                    // methods
            .def("assign_options",[] (CommandLineParser  & ptr, CommandLineOptions const & cl_options, Options & options, String prefix) {
                ptr.assign_options(cl_options, options, prefix); } )
            ;

    py::class_<EnvManager, std::shared_ptr<EnvManager>>(m, "EnvManager")
            // ctors
            .def(py::init<Strings const &>())
                    // methods
            .def("add_env",[] (EnvManager  & ptr, String const & env) {
                ptr.add_env(env); } )
            .def("set_envs",[] (EnvManager  & ptr) {
                ptr.set_envs(); } )
            ;



    py::class_<Option, std::shared_ptr<Option>>(m, "Option")
            // ctors
            .def(py::init<>())
            .def(py::init<Option const &>())
            .def(py::init<String const &,float const &,OptionType const &>())
            .def(py::init<String const &,String const &,OptionType const &>())
            .def(py::init<String const &,char const *,OptionType const &>())
            .def(py::init<char const *,char const *,OptionType const &>())
            .def(py::init<String const &,bool const &,OptionType const &>())
            .def(py::init<String const &,int const &,OptionType const &>())
                    // methods
            .def("get_float",[] (Option  & ptr) -> float {
                return ptr.get_float(); } )
            .def("get_int",[] (Option  & ptr) -> int {
                return ptr.get_int(); } )
            .def("get_string",[] (Option  & ptr) -> String const & {
                return ptr.get_string(); } )
            .def("get_bool",[] (Option  & ptr) -> bool const & {
                return ptr.get_bool(); } )
            .def("value",[] (Option  & ptr, int const & i) {
                ptr.value(i); } )
            .def("value",[] (Option  & ptr, char const * c_str) {
                ptr.value(c_str); } )
            .def("value",[] (Option  & ptr, int const & i) {
                ptr.value(i); } )
            .def("value",[] (Option  & ptr, int const & i) {
                ptr.value(i); } )
            .def("value",[] (Option  & ptr, int const & i) {
                ptr.value(i); } )
            .def("name",[] (Option const & ptr) ->  String const & {
                return ptr.name(); } )
            .def("type",[] (Option const & ptr) ->  OptionType const & {
                return ptr.type(); } )
            .def("type_name",[] (Option  & ptr) ->  String {
                return ptr.type_name(); } )
            ;

    py::class_<OptionClass, std::shared_ptr<OptionClass>>(m, "OptionClass")
            ;

    py::class_<OptionType, std::shared_ptr<OptionType>>(m, "OptionType")
            ;

    py::class_<Options, std::shared_ptr<Options>>(m, "Options")
            // ctors
            .def(py::init<>())
            .def(py::init<String const &>())
            .def(py::init<Options const &>())
                    // methods
            .def("begin",[] (Options  & ptr) -> std::vector<OptionOP>::iterator {
                return ptr.begin(); } )
            .def("end",[] (Options  & ptr) -> std::vector<OptionOP>::iterator {
                return ptr.end(); } )
            .def("begin",[] (Options const & ptr) -> std::vector<OptionOP>::const_iterator {
                return ptr.begin(); } )
            .def("end",[] (Options const & ptr) -> std::vector<OptionOP>::const_iterator {
                return ptr.end(); } )
            .def("size",[] (Options  & ptr) -> size_t {
                return ptr.size(); } )
            .def("lock_option_adding",[] (Options  & ptr) {
                ptr.lock_option_adding(); } )
            .def("add_option",[] (Options  & ptr, String const & name, int const & val, OptionType const & type) { ptr.add_option(name, val, type); } )
            .def("add_option",[] (Options  & ptr, String const & name, bool const & val, OptionType const & type) { ptr.add_option(name, val, type); } )
            .def("add_option",[] (Options  & ptr, String const & name, String const & val, OptionType const & type) { ptr.add_option(name, val, type); } )
            .def("add_option",[] (Options  & ptr, String const & name, float const & val, OptionType const & type) { ptr.add_option(name, val, type); } )
            .def("get_int",[] (Options  & ptr, String const & name) ->  int {
                return ptr.get_int(name); } )
            .def("get_float",[] (Options  & ptr, String const & name) ->  float {
                return ptr.get_float(name); } )
            .def("get_string",[] (Options  & ptr, String const & name) ->  String const & {
                return ptr.get_string(name); } )
            .def("get_bool",[] (Options  & ptr, String const & name) ->  bool {
                return ptr.get_bool(name); } )
            .def("has_option",[] (Options  & ptr, String const & name) ->  bool {
                return ptr.has_option(name); } )
            .def("set_value",[] (Options  & ptr, String const & name, bool const & val) { ptr.set_value(name, val); } )
            .def("set_value",[] (Options  & ptr, String const & name, float const & val) { ptr.set_value(name, val); } )
            .def("set_value",[] (Options  & ptr, String const & name, int const & val) { ptr.set_value(name, val); } )
            .def("set_value",[] (Options  & ptr, String const & name, String const & val) { ptr.set_value(name, val); } )
            ;

    py::class_<RNAMakeIOException, std::shared_ptr<RNAMakeIOException>>(m, "RNAMakeIOException")
            // ctors
            .def(py::init<String const &>())
            .def(py::init<const char *>())
            ;

    py::class_<RNAMakeImplementationExcepetion, std::shared_ptr<RNAMakeImplementationExcepetion>>(m, "RNAMakeImplementationExcepetion")
            // ctors
            .def(py::init<String const &>())
            .def(py::init<const char *>())
            ;
/*
        py::class_<VectorContainer, std::shared_ptr<VectorContainer>>(m, "VectorContainer")
		// ctors
		.def(py::init<std::vector<T> const &>())
		// methods
		.def("size",[] (VectorContainer  & ptr) ->  size_t {
		 return ptr.size(); } )
		.def("begin",[] (VectorContainer const & ptr) -> const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (VectorContainer const & ptr) -> const_iterator {
		 return ptr.end(); } )
		.def("at",[] (VectorContainer const & ptr, Index i) ->  T const & {
		 return ptr.at(i); } )
		.def("get_data",[] (VectorContainer  & ptr) ->  std::vector<T> const & {
		 return ptr.get_data(); } )
		// operators
		.def(py::self [] Index)
		;
*/
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// math
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // free functions

    m.def("are_floats_equal", [] (double const a, double const b, double tol) -> int {
              return are_floats_equal(a, b, tol); },
          py::arg("a"),
          py::arg("b"),
          py::arg("tol") = 0.001
    );

    m.def("are_xyzMatrix_equal", [] (Matrix const & m, Matrix const & mc) -> int {
              return are_xyzMatrix_equal(m, mc); },
          py::arg("m"),
          py::arg("mc")
    );

    m.def("are_xyzVector_equal", [] (Vector const & vec, Vector const & correct_vec, float tol) -> int {
              return are_xyzVector_equal(vec, correct_vec, tol); },
          py::arg("vec"),
          py::arg("correct_vec"),
          py::arg("tol") = 0.001
    );

    m.def("are_xyzVectors_equal", [] (Vectors const & v, Vectors const & vc) -> int {
              return are_xyzVectors_equal(v, vc); },
          py::arg("v"),
          py::arg("vc")
    );

    m.def("avg_unsigned_diff", [] (std::vector<double> const & x, std::vector<double> const & y) -> double {
              return avg_unsigned_diff(x, y); },
          py::arg("x"),
          py::arg("y")
    );

    m.def("axis_angle_from_matrix", [] (Matrix & m, AxisAngle & aa) {
              axis_angle_from_matrix(m, aa); },
          py::arg("m"),
          py::arg("aa")
    );

    m.def("calc_euler", [] (Matrix & M, Vector & euler) {
              calc_euler(M, euler); },
          py::arg("M"),
          py::arg("euler")
    );

    m.def("degrees", [] (float radians) -> float {
              return degrees(radians); },
          py::arg("radians")
    );

    m.def("dot_vector", [] (math::Matrix m, Vector const & v, Vector & vr) {
              dot_vector(m, v, vr); },
          py::arg("m"),
          py::arg("v"),
          py::arg("vr")
    );

    m.def("dot_vectors", [] (math::Matrix m, Vectors const & v, Vectors & vr) {
              dot_vectors(m, v, vr); },
          py::arg("m"),
          py::arg("v"),
          py::arg("vr")
    );

    m.def("get_quaternion_from_matrix", [] (Matrix const & m) -> Quaternion {
              return get_quaternion_from_matrix(m); },
          py::arg("m")
    );

    m.def("get_random_quaternion", [] () -> Quaternion {
        return get_random_quaternion(); }
    );

    m.def("matrix_from_str", [] (std::string const & s) ->  const Matrix {
              return matrix_from_str(s); },
          py::arg("s")
    );

    m.def("matrix_to_str", [] (math::Matrix m) ->  String {
              return matrix_to_str(m); },
          py::arg("m")
    );

    m.def("mean", [] (std::vector<double> const & a) -> double {
              return mean(a); },
          py::arg("a")
    );

    m.def("norm", [] (std::vector<double> const & v) ->  double {
              return norm(v); },
          py::arg("v")
    );

    m.def("pearson_coeff", [] (std::vector<double> const & x, std::vector<double> const & y) -> double {
              return pearson_coeff(x, y); },
          py::arg("x"),
          py::arg("y")
    );

    m.def("power_iteration", [] (std::vector<std::vector<double> > const & A, std::vector<double> & eigen_values, int num_simulations) {
              power_iteration(A, eigen_values, num_simulations); },
          py::arg("A"),
          py::arg("eigen_values"),
          py::arg("num_simulations")
    );

    m.def("sqsum", [] (std::vector<double> const & a) -> double {
              return sqsum(a); },
          py::arg("a")
    );

    m.def("stdev", [] (std::vector<double> const & nums) -> double {
              return stdev(nums); },
          py::arg("nums")
    );

    m.def("sum", [] (std::vector<double> const & a) -> double {
              return sum(a); },
          py::arg("a")
    );

    m.def("vector_from_str", [] (std::string const & s) ->  const Vector {
              return vector_from_str(s); },
          py::arg("s")
    );

    m.def("vector_to_str", [] (math::Vector v) ->  const String {
              return vector_to_str(v); },
          py::arg("v")
    );

    m.def("vectors_from_str", [] (std::string const & s) ->  const Vectors {
              return vectors_from_str(s); },
          py::arg("s")
    );

    m.def("vectors_to_str", [] (math::Vectors vs) ->  String {
              return vectors_to_str(vs); },
          py::arg("vs")
    );
    // classes

    py::class_<AverageQuaternionCalculator, std::shared_ptr<AverageQuaternionCalculator>>(m, "AverageQuaternionCalculator")
            // ctors
            .def(py::init<>())
                    // methods
            .def("add_quaternion",[] (AverageQuaternionCalculator  & ptr, Quaternion const & q) {
                ptr.add_quaternion(q); } )
            .def("get_average",[] (AverageQuaternionCalculator  & ptr) -> Quaternion {
                return ptr.get_average(); } )
            ;

    py::class_<AxisAngle, std::shared_ptr<AxisAngle>>(m, "AxisAngle")
            // public attributes
            .def_readwrite("angle", &AxisAngle::angle)
            .def_readwrite("axis", &AxisAngle::axis)
            ;

    py::class_<Quaternion, std::shared_ptr<Quaternion>>(m, "Quaternion")
            // ctors
            .def(py::init<>())
            .def(py::init<double>())
            .def(py::init<double,double,double,double>())
            .def(py::init<Quaternion const &>())
                    // methods
            .def("dot",[] (Quaternion const & ptr, Quaternion const & q) ->  double {
                return ptr.dot(q); } )
            .def("get_rotation_matrix",[] (Quaternion  & ptr) -> Matrix {
                return ptr.get_rotation_matrix(); } )
            .def("get_a",[] (Quaternion const & ptr) ->  double {
                return ptr.get_a(); } )
            .def("get_b",[] (Quaternion const & ptr) ->  double {
                return ptr.get_b(); } )
            .def("get_c",[] (Quaternion const & ptr) ->  double {
                return ptr.get_c(); } )
            .def("get_d",[] (Quaternion const & ptr) ->  double {
                return ptr.get_d(); } )
            .def("to_str", [] (Quaternion const& q) -> String { auto ss = std::stringstream(); ss << q; return ss.str(); })
            .def("at", [](Quaternion const& q, int index) -> decltype(q[index]) {return q[index];})
                    // operators
                    //.def(py::self << std::ostream)
                    //.def(py::self [] int())
            .def(py::self += double())
            .def(py::self *= double())
            ;

    py::class_<SixDCoordinateBinner, std::shared_ptr<SixDCoordinateBinner>>(m, "SixDCoordinateBinner")
            // ctors
            .def(py::init<math::BoundingBox,math::Real6>())
                    // methods
            .def("bin6",[] (SixDCoordinateBinner const & ptr, math::Real6 values) -> Bin6D {
                return ptr.bin6(values); } )
            .def("bin_center_point",[] (SixDCoordinateBinner const & ptr, math::Bin6D bin) -> Real6 {
                return ptr.bin_center_point(bin); } )
            .def("bin_to_values",[] (SixDCoordinateBinner const & ptr, math::Bin6D bin) -> Real6 {
                return ptr.bin_to_values(bin); } )
            .def("bin_index",[] (SixDCoordinateBinner const & ptr, math::Real6 values) -> uint64_t {
                return ptr.bin_index(values); } )
            .def("bin_from_index",[] (SixDCoordinateBinner  & ptr, uint64_t bin_index) -> Bin6D {
                return ptr.bin_from_index(bin_index); } )
            .def("get_bounding_box",[] (SixDCoordinateBinner  & ptr) -> BoundingBox const & {
                return ptr.get_bounding_box(); } )
            .def("get_bin_widths",[] (SixDCoordinateBinner  & ptr) -> Real6 const & {
                return ptr.get_bin_widths(); } )
            ;

    py::class_<SixDHistogram, std::shared_ptr<SixDHistogram>>(m, "SixDHistogram")
            // ctors
            .def(py::init<math::BoundingBox,math::Real6>())
            .def(py::init<Strings const &,SixDHistogramStrType const &>())
            .def(py::init<std::ifstream &>())
                    // methods
                    //.def("read",[] (SixDHistogram  & ptr ) { return ptr.read(>(&num), ); } )
                    //.def("read",[] (SixDHistogram  & ptr, reinterpret_cast<char * >(&count), sizeof(count) ) -> in {
                    //return ptr.read(>(&count), ); } )
            .def("size",[] (SixDHistogram  & ptr) ->  size_t {
                return ptr.size(); } )
            .def("add",[] (SixDHistogram  & ptr, math::Real6 values) {
                ptr.add(values); } )
            .def("contains",[] (SixDHistogram  & ptr, math::Real6 values) -> bool {
                return ptr.contains(values); } )
            .def("get",[] (SixDHistogram  & ptr, math::Real6 values) -> uint64_t {
                return ptr.get(values); } )
            .def("within_constraints",[] (SixDHistogram  & ptr, std::array<Real2, 6> const & constraints) -> uint64_t {
                return ptr.within_constraints(constraints); } )
            .def("total_count",[] (SixDHistogram  & ptr) -> uint64_t {
                return ptr.total_count(); } )
            .def("to_text_file",[] (SixDHistogram  & ptr, String const & fname) {
                ptr.to_text_file(fname); } )
            .def("to_binary_file",[] (SixDHistogram  & ptr, String const & fname, uint64_t cuttoff) {
                ptr.to_binary_file(fname, cuttoff); } )
            .def("output_binary",[] (SixDHistogram  & ptr, std::ofstream & out, uint64_t cuttoff) {
                ptr.output_binary(out, cuttoff); } )
        // public attributes
        //.def_readwrite("datas", &SixDHistogram::)
        ////.def_readwrite("", &SixDHistogram::)
        //.def_readwrite("lower", &SixDHistogram::lower)
        //.def_readwrite("", &SixDHistogram::)
        //.def_readwrite("upper", &SixDHistogram::upper)
        //.def_readwrite("bin_widths", &SixDHistogram::bin_widths)
        //.def_readwrite("", &SixDHistogram::)
        //.def_readwrite("binner_", &SixDHistogram::binner_)
        //.def_readwrite("num", &sum)
        //.def_readwrite("key", &SixDHistogram::key)
        //.def_readwrite("count", &SixDHistogram::count)
        //.def_readwrite("num", &SixDHistogram::num)
        //.def_readwrite("", &SixDHistogram::)
            ;
//        py::class_<SixDHistogramStrType, std::shared_ptr<SixDHistogramStrType>>(m, "SixDHistogramStrType")
//		;
    py::enum_<SixDHistogramStrType>(m, "SixDHistogramStrType")
            .value("TEXT", SixDHistogramStrType::TEXT)
            .value("BINARY", SixDHistogramStrType::BINARY)
            ;

    py::class_<ThreeDCoordinateBinner, std::shared_ptr<ThreeDCoordinateBinner>>(m, "ThreeDCoordinateBinner")
            // ctors
            .def(py::init<math::BoundingBox,math::Real3>())
                    // methods
            .def("bin3",[] (ThreeDCoordinateBinner const & ptr, Point const & values) -> Real3 {
                return ptr.bin3(values); } )
            .def("bin_center_point",[] (ThreeDCoordinateBinner const & ptr, math::Real3 bin) -> Real3 {
                return ptr.bin_center_point(bin); } )
            .def("bin_to_values",[] (ThreeDCoordinateBinner const & ptr, math::Bin6D bin) -> Real3 {
                return ptr.bin_to_values(bin); } )
            .def("bin_index",[] (ThreeDCoordinateBinner const & ptr, Point const & values) -> uint32_t {
                return ptr.bin_index(values); } )
            .def("bin_from_index",[] (ThreeDCoordinateBinner  & ptr, uint32_t bin_index) -> Real3 {
                return ptr.bin_from_index(bin_index); } )
            .def("get_bounding_box",[] (ThreeDCoordinateBinner  & ptr) -> BoundingBox const & {
                return ptr.get_bounding_box(); } )
            .def("get_bin_widths",[] (ThreeDCoordinateBinner  & ptr) -> Real3 const & {
                return ptr.get_bin_widths(); } )
            ;

    py::class_<ThreeDHistogram, std::shared_ptr<ThreeDHistogram>>(m, "ThreeDHistogram")
            // ctors
            .def(py::init<>())
            .def(py::init<math::BoundingBox,math::Real3>())
                    // methods
            .def("setup",[] (ThreeDHistogram  & ptr, math::BoundingBox bounding_box, math::Real3 bin_widths) {
                ptr.setup(bounding_box, bin_widths); } )
            .def("size",[] (ThreeDHistogram  & ptr) ->  size_t {
                return ptr.size(); } )
            .def("add",[] (ThreeDHistogram  & ptr, Point const & values) {
                ptr.add(values); } )
            .def("contains",[] (ThreeDHistogram  & ptr, Point const & values) -> bool {
                return ptr.contains(values); } )
            .def("get",[] (ThreeDHistogram  & ptr, Point const & values) -> uint32_t {
                return ptr.get(values); } )
            .def("write_histo_to_pdb",[] (ThreeDHistogram  & ptr, String const & pdb_name) {
                ptr.write_histo_to_pdb(pdb_name); } )
            ;

    py::class_<Transform, std::shared_ptr<Transform>>(m, "Transform")
            // ctors
            .def(py::init<>())
            .def(py::init<Matrix const &,Point const &>())
                    // methods
                    //.def("dot",[] (Transform  & ptr, Matrix const & a, Matrix const & b, Transform & c) { ptr.dot(a, b, c); } )
            .def("xx",[] (Transform const & ptr) -> float {
                return ptr.xx(); } )
            .def("xy",[] (Transform const & ptr) -> float {
                return ptr.xy(); } )
            .def("xz",[] (Transform const & ptr) -> float {
                return ptr.xz(); } )
            .def("yx",[] (Transform const & ptr) -> float {
                return ptr.yx(); } )
            .def("yy",[] (Transform const & ptr) -> float {
                return ptr.yy(); } )
            .def("yz",[] (Transform const & ptr) -> float {
                return ptr.yz(); } )
            .def("zx",[] (Transform const & ptr) -> float {
                return ptr.zx(); } )
            .def("zy",[] (Transform const & ptr) -> float {
                return ptr.zy(); } )
            .def("zz",[] (Transform const & ptr) -> float {
                return ptr.zz(); } )
            .def("px",[] (Transform const & ptr) -> float {
                return ptr.px(); } )
            .def("py",[] (Transform const & ptr) -> float {
                return ptr.py(); } )
            .def("pz",[] (Transform const & ptr) -> float {
                return ptr.pz(); } )
            .def("xx",[] (Transform  & ptr) -> float {
                return ptr.xx(); } )
            .def("xy",[] (Transform  & ptr) -> float {
                return ptr.xy(); } )
            .def("xz",[] (Transform  & ptr) -> float {
                return ptr.xz(); } )
            .def("yx",[] (Transform  & ptr) -> float {
                return ptr.yx(); } )
            .def("yy",[] (Transform  & ptr) -> float {
                return ptr.yy(); } )
            .def("yz",[] (Transform  & ptr) -> float {
                return ptr.yz(); } )
            .def("zx",[] (Transform  & ptr) -> float {
                return ptr.zx(); } )
            .def("zy",[] (Transform  & ptr) -> float {
                return ptr.zy(); } )
            .def("zz",[] (Transform  & ptr) -> float {
                return ptr.zz(); } )
            .def("px",[] (Transform  & ptr) -> float {
                return ptr.px(); } )
            .def("py",[] (Transform  & ptr) -> float {
                return ptr.py(); } )
            .def("pz",[] (Transform  & ptr) -> float {
                return ptr.pz(); } )
            .def("xaxis",[] (Transform const & ptr) -> Vector {
                return ptr.xaxis(); } )
            .def("yaxis",[] (Transform const & ptr) -> Vector {
                return ptr.yaxis(); } )
            .def("zaxis",[] (Transform const & ptr) -> Vector {
                return ptr.zaxis(); } )
            .def("translation",[] (Transform const & ptr) -> Vector {
                return ptr.translation(); } )
            .def("rotation",[] (Transform const & ptr) ->  Matrix {
                return ptr.rotation(); } )
            .def("rotation",[] (Transform  & ptr, Matrix const & m) {
                ptr.rotation(m); } )
            .def("translation",[] (Transform  & ptr, Vector const & v) {
                ptr.translation(v); } )
            ;

    py::class_<_BoundingBox<Point>, std::shared_ptr<_BoundingBox<Point>>>(m, "_BoundingBox")
            // ctors
            .def(py::init<>())
            .def(py::init<Point const &>())
            .def(py::init<Point const &, Point const &>())
            .def(py::init<_BoundingBox<Point> const &>())
                    // methods
//		.def("add",[] (_BoundingBox<Point>  & ptr, Point const & pp) {
//		ptr.add(pp); } )
            .def("reset",[] (_BoundingBox<Point>  & ptr, Point const & p) {
                ptr.reset(p); } )
            .def("expand",[] (_BoundingBox<Point>  & ptr, double const & scalar) {
                ptr.expand(scalar); } )
            .def("contract",[] (_BoundingBox<Point>  & ptr, double const & scalar) {
                ptr.contract(scalar); } )
            .def("translate",[] (_BoundingBox<Point>  & ptr, Point const & t) {
                ptr.translate(t); } )
            .def("intersects",[] (_BoundingBox<Point> const & ptr, _BoundingBox<Point> const & bb) ->  bool {
                return ptr.intersects(bb); } )
            .def("contains",[] (_BoundingBox<Point> const & ptr, double const & x, double const & y, double const & z) ->  bool {
                return ptr.contains(x, y, z); } )
            .def("contains",[] (_BoundingBox<Point> const & ptr, Point const & p) ->  bool {
                return ptr.contains(p); } )
            .def("set_lower",[] (_BoundingBox<Point>  & ptr, Point const & p) {
                ptr.set_lower(p); } )
            .def("set_upper",[] (_BoundingBox<Point>  & ptr, Point const & p) {
                ptr.set_upper(p); } )
            .def("lower",[] (_BoundingBox<Point> const & ptr) ->  Point const & {
                return ptr.lower(); } )
            .def("upper",[] (_BoundingBox<Point> const & ptr) ->  Point const & {
                return ptr.upper(); } )
        // operators
        //.def(py::self = py::self)
            ;

#define Value double
    py::class_<xyzMatrix<Value>, std::shared_ptr<xyzMatrix<Value>>>(m, "xyzMatrix")
            // ctors
            .def(py::init<>())
            .def(py::init<xyzMatrix<Value> const &>())
            .def(py::init<xyzMatrix<double> const &>())
            .def(py::init<const double &,const double &,const double &,const double &,const double &,const double &,const double &,const double &,const double &>())
            .def(py::init<Value const &>())
            .def(py::init<String const &>())
                    // methods
            .def("transpose",[] (xyzMatrix<Value> & ptr) { return ptr.transpose(); } )
                    //.def("dot",[] (xyzMatrix<Value>  & ptr, xyzMatrix<Value> const & a, xyzMatrix<Value> const & b, xyzMatrix<Value> & c) {
                    //ptr.dot(a, b, c); } )
            .def("to_str",[] (xyzMatrix<Value> const & ptr) ->  String const {
                return ptr.to_str(); } )
            .def("row",[] (xyzMatrix<Value>  & ptr, const int i, Vector const & v) {
                ptr.row(i, v); } )
            .def("row",[] (xyzMatrix<Value>  & ptr, const int i, std::vector<Value> const & v) {
                ptr.row(i, v); } )
            .def("identity",[] (xyzMatrix<Value>  & ptr) ->  xyzMatrix<Value> {
                return ptr.identity(); } )
            .def("transpose",[] (xyzMatrix<Value>  & ptr) ->  xyzMatrix<Value> & {
                return ptr.transpose(); } )
            .def("difference",[] (xyzMatrix<Value> const & ptr, xyzMatrix<Value> const & b) ->  float const {
                return ptr.difference(b); } )
            .def("get_flip_orientation",[] (xyzMatrix<Value> const & ptr) ->  const xyzMatrix<Value> {
                return ptr.get_flip_orientation(); } )
            .def("get_unitarize",[] (xyzMatrix<Value> const & ptr) ->  const xyzMatrix<Value> {
                return ptr.get_unitarize(); } )
            .def("unitarize",[] (xyzMatrix<Value>  & ptr) {
                ptr.unitarize(); } )
            .def("xx",[] (xyzMatrix<Value> const & ptr) ->  Value const & {
                return ptr.xx(); } )
            .def("xy",[] (xyzMatrix<Value> const & ptr) ->  Value const & {
                return ptr.xy(); } )
            .def("xz",[] (xyzMatrix<Value> const & ptr) ->  Value const & {
                return ptr.xz(); } )
            .def("yx",[] (xyzMatrix<Value> const & ptr) ->  Value const & {
                return ptr.yx(); } )
            .def("yy",[] (xyzMatrix<Value> const & ptr) ->  Value const & {
                return ptr.yy(); } )
            .def("yz",[] (xyzMatrix<Value> const & ptr) ->  Value const & {
                return ptr.yz(); } )
            .def("zx",[] (xyzMatrix<Value> const & ptr) ->  Value const & {
                return ptr.zx(); } )
            .def("zy",[] (xyzMatrix<Value> const & ptr) ->  Value const & {
                return ptr.zy(); } )
            .def("zz",[] (xyzMatrix<Value> const & ptr) ->  Value const & {
                return ptr.zz(); } )
            .def("xx",[] (xyzMatrix<Value>  & ptr, Value const & xx_a) {
                ptr.xx(xx_a); } )
            .def("xy",[] (xyzMatrix<Value>  & ptr, Value const & xy_a) {
                ptr.xy(xy_a); } )
            .def("xz",[] (xyzMatrix<Value>  & ptr, Value const & xz_a) {
                ptr.xz(xz_a); } )
            .def("yx",[] (xyzMatrix<Value>  & ptr, Value const & yx_a) {
                ptr.yx(yx_a); } )
            .def("yy",[] (xyzMatrix<Value>  & ptr, Value const & yy_a) {
                ptr.yy(yy_a); } )
            .def("yz",[] (xyzMatrix<Value>  & ptr, Value const & yz_a) {
                ptr.yz(yz_a); } )
            .def("zx",[] (xyzMatrix<Value>  & ptr, Value const & zx_a) {
                ptr.zx(zx_a); } )
            .def("zy",[] (xyzMatrix<Value>  & ptr, Value const & zy_a) {
                ptr.zy(zy_a); } )
            .def("zz",[] (xyzMatrix<Value>  & ptr, Value const & zz_a) {
                ptr.zz(zz_a); } )
            .def("transposed",[] (xyzMatrix<Value> const & ptr) ->  xyzMatrix<Value> {
                return ptr.transposed(); } )
                    // operators
            .def(py::self += py::self )
            .def(py::self -= py::self )
            .def(py::self += Value() )
            .def(py::self -= Value() )
            .def(py::self + py::self)
            .def(py::self + py::self)
            .def(py::self + Value() )
            .def(py::self - py::self)
            .def(py::self - py::self)
            .def(py::self - Value())
            .def(py::self * py::self)
            .def(py::self * py::self)
            .def(py::self * Value())
        //.def(py::self / py::self)
            ;

    py::class_<xyzVector<Value>, std::shared_ptr<xyzVector<Value>>>(m, "xyzVector<Value>")
            // ctors
            .def(py::init<>())
            .def(py::init<xyzVector<Value> const &>())
            .def(py::init<xyzVector<Value> const &>())
            .def(py::init<Value const &,Value const &,Value const &>())
            .def(py::init<std::vector<Value> const &>())
            .def(py::init<Value const &>())
            .def(py::init<String const &>())
                    // methods
            .def("to_str",[] (xyzVector<Value> const & ptr) ->  String const {
                return ptr.to_str(); } )
            .def("zero",[] (xyzVector<Value>  & ptr) ->  xyzVector<Value> & {
                return ptr.zero(); } )
            .def("negate",[] (xyzVector<Value>  & ptr) ->  xyzVector<Value> & {
                return ptr.negate(); } )
            .def("negated",[] (xyzVector<Value> const & ptr) ->  xyzVector<Value> {
                return ptr.negated(); } )
            .def("negated",[] (xyzVector<Value> const & ptr, xyzVector<Value> & a) {
                ptr.negated(a); } )
//		.def("add",[] (xyzVector<Value>  & ptr, xyzVector<Value> const & a, xyzVector<Value> const & b, xyzVector<Value> & r) {
//		ptr.add(a, b, r); } )
//		.def("add",[] (xyzVector<Value>  & ptr, xyzVector<Value> const & v, Value const & t, xyzVector<Value> & r) {
//		ptr.add(v, t, r); } )
//		.def("add",[] (xyzVector<Value>  & ptr, Value const & t, xyzVector<Value> const & v, xyzVector<Value> & r) {
//		ptr.add(t, v, r); } )
//		.def("subtract",[] (xyzVector<Value>  & ptr, xyzVector<Value> const & a, xyzVector<Value> const & b, xyzVector<Value> & r) {
//		ptr.subtract(a, b, r); } )
//		.def("subtract",[] (xyzVector<Value>  & ptr, xyzVector<Value> const & v, Value const & t, xyzVector<Value> & r) {
//		ptr.subtract(v, t, r); } )
//		.def("subtract",[] (xyzVector<Value>  & ptr, Value const & t, xyzVector<Value> const & v, xyzVector<Value> & r) {
//		ptr.subtract(t, v, r); } )
//		.def("multiply",[] (xyzVector<Value>  & ptr, xyzVector<Value> const & v, Value const & t, xyzVector<Value> & r) {
//		ptr.multiply(v, t, r); } )
//		.def("multiply",[] (xyzVector<Value>  & ptr, Value const & t, xyzVector<Value> const & v, xyzVector<Value> & r) {
//		ptr.multiply(t, v, r); } )
//		.def("divide",[] (xyzVector<Value>  & ptr, xyzVector<Value> const & v, Value const & t, xyzVector<Value> & r) {
//		ptr.divide(v, t, r); } )
            .def("normalize",[] (xyzVector<Value>  & ptr) ->  xyzVector<Value> & {
                return ptr.normalize(); } )
            .def("distance",[] (xyzVector<Value> const & ptr, xyzVector<Value> const & v) ->  Value {
                return ptr.distance(v); } )
            .def("distance_squared",[] (xyzVector<Value> const & ptr, xyzVector<Value> const & v) ->  Value {
                return ptr.distance_squared(v); } )
            .def("dot",[] (xyzVector<Value> const & ptr, xyzVector<Value> const & v) ->  Value {
                return ptr.dot(v); } )
//		.def("dot_product",[] (xyzVector<Value>  & ptr, xyzVector<Value> const & a, xyzVector<Value> const & b) -> friend  Value {
//		 return ptr.dot_product(a, b); } )
            .def("cross",[] (xyzVector<Value> const & ptr, xyzVector<Value> const & v) ->  xyzVector<Value> {
                return ptr.cross(v); } )
//		.def("cross",[] (xyzVector<Value>  & ptr, xyzVector<Value> const & a, xyzVector<Value> const & b) -> friend  xyzVector<Value> {
//		 return ptr.cross(a, b); } )
            .def("x",[] (xyzVector<Value> const & ptr) ->  Value const & {
                return ptr.x(); } )
            .def("y",[] (xyzVector<Value> const & ptr) ->  Value const & {
                return ptr.y(); } )
            .def("z",[] (xyzVector<Value> const & ptr) ->  Value const & {
                return ptr.z(); } )
            .def("length",[] (xyzVector<Value> const & ptr) ->  Value {
                return ptr.length(); } )
            .def("length_squared",[] (xyzVector<Value> const & ptr) ->  Value {
                return ptr.length_squared(); } )
            .def("norm",[] (xyzVector<Value> const & ptr) ->  Value {
                return ptr.norm(); } )
            .def("norm_squared",[] (xyzVector<Value> const & ptr) ->  Value {
                return ptr.norm_squared(); } )
            .def("magnitude",[] (xyzVector<Value> const & ptr) ->  Value {
                return ptr.magnitude(); } )
            .def("magnitude_squared",[] (xyzVector<Value> const & ptr) ->  Value {
                return ptr.magnitude_squared(); } )
            .def("x",[] (xyzVector<Value>  & ptr, Value const & x_a) {
                ptr.x(x_a); } )
            .def("y",[] (xyzVector<Value>  & ptr, Value const & y_a) {
                ptr.y(y_a); } )
            .def("z",[] (xyzVector<Value>  & ptr, Value const & z_a) {
                ptr.z(z_a); } )
                    // operators
//		.def(py::self = py::self)
//		.def(py::self = py::self  )
            .def(py::self += py::self  )
            .def(py::self -= py::self  )
                    //.def(py::self = Value()  )
            .def(py::self += Value()  )
            .def(py::self -= Value()  )
            .def(py::self *= Value()  )
            .def(py::self /= Value()  )
            .def(py::self + py::self)
            .def(py::self + py::self)
            .def(py::self + Value()  )
            .def(py::self - py::self)
            .def(py::self - py::self)
            .def(py::self - Value()  )
//		.def(py::self * py::self)
            .def(py::self * Value()  )
//		.def(py::self / py::self)
            .def("at", [] (xyzVector<Value> const& vec, int index) -> Value { return vec[index]; })
//        .def(py::self [] int )
//		.def(py::self [] int )
//		.def(py::self () int )
//		.def(py::self () int )
            .def(py::self == py::self)
            .def(py::self != py::self)
            ;
#undef Value
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// data_structure
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
        py::class_<AdjacencyList, std::shared_ptr<AdjacencyList>>(m, "AdjacencyList")
		// ctors
		.def(py::init<>())
		.def(py::init<AdjacencyList const &>())
		// methods
		.def("begin",[] (AdjacencyList const & ptr) -> const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (AdjacencyList const & ptr) -> const_iterator {
		 return ptr.end(); } )
		.def("add_node",[] (AdjacencyList  & ptr, DataType const & d, Size n_edges) -> Index {
		 return ptr.add_node(d, n_edges); } )
		.def("add_edge",[] (AdjacencyList  & ptr, NodeIndexandEdge const & nie1, NodeIndexandEdge const & nie2) {
		ptr.add_edge(nie1, nie2); } )
		.def("remove_node",[] (AdjacencyList  & ptr, Index ni) {
		ptr.remove_node(ni); } )
		.def("remove_edge",[] (AdjacencyList  & ptr, NodeIndexandEdge const & nie1, NodeIndexandEdge const & nie2) {
		ptr.remove_edge(nie1, nie2); } )
		.def("get_num_nodes",[] (AdjacencyList const & ptr) ->  size_t {
		 return ptr.get_num_nodes(); } )
		.def("get_num_edges",[] (AdjacencyList const & ptr) -> size_t {
		 return ptr.get_num_edges(); } )
		.def("get_node_edges",[] (AdjacencyList const & ptr, Index ni) -> std::vector<Edge const *> const & {
		 return ptr.get_node_edges(ni); } )
		.def("get_node_data",[] (AdjacencyList const & ptr, Index ni) -> DataType const & {
		 return ptr.get_node_data(ni); } )
		.def("get_node_data",[] (AdjacencyList  & ptr, Index ni) -> DataType & {
		 return ptr.get_node_data(ni); } )
		.def("get_node",[] (AdjacencyList const & ptr, Index ni) -> Node<DataType> const & {
		 return ptr.get_node(ni); } )
		.def("get_node",[] (AdjacencyList  & ptr, Index ni) -> Node<DataType> & {
		 return ptr.get_node(ni); } )
		.def("get_connected_node_info",[] (AdjacencyList const & ptr, NodeIndexandEdge const & nei) -> NodeIndexandEdge {
		 return ptr.get_connected_node_info(nei); } )
		.def("edge_between_nodes",[] (AdjacencyList const & ptr, Index n1, Index n2) -> bool {
		 return ptr.edge_between_nodes(n1, n2); } )
		.def("edge_index_empty",[] (AdjacencyList const & ptr, Index ni, Index ei) -> bool {
		 return ptr.edge_index_empty(ni, ei); } )
		// operators
		.def(py::self = py::self)
		;

        py::class_<DirectedAdjacencyList, std::shared_ptr<DirectedAdjacencyList>>(m, "DirectedAdjacencyList")
		// ctors
		.def(py::init<>())
		.def(py::init<DirectedAdjacencyList const &>())
		// methods
		.def("add_node",[] (DirectedAdjacencyList  & ptr, DataType const & d, Size n_edges) -> Index {
		 return ptr.add_node(d, n_edges); } )
		.def("add_node",[] (DirectedAdjacencyList  & ptr, DataType const & d, Size n_edges, Index n_end_index, NodeIndexandEdge const & pie) -> Index {
		 return ptr.add_node(d, n_edges, n_end_index, pie); } )
		.def("remove_node",[] (DirectedAdjacencyList  & ptr, Index ni) {
		ptr.remove_node(ni); } )
		.def("has_parent",[] (DirectedAdjacencyList const & ptr, Index n_index) -> bool {
		 return ptr.has_parent(n_index); } )
		.def("get_parent_index",[] (DirectedAdjacencyList const & ptr, Index n_index) -> Index {
		 return ptr.get_parent_index(n_index); } )
		.def("get_parent_end_index",[] (DirectedAdjacencyList const & ptr, Index n_index) -> Index {
		 return ptr.get_parent_end_index(n_index); } )
		// operators
		.def(py::self = py::self)
		;

        py::class_<DirectedGraph, std::shared_ptr<DirectedGraph>>(m, "DirectedGraph")
		// ctors
		.def(py::init<>())
		.def(py::init<DirectedGraph const &>())
		// methods
		.def("setup_sub_graph_transversal",[] (DirectedGraph  & ptr, Index start_n, Index end_n) {
		ptr.setup_sub_graph_transversal(start_n, end_n); } )
		.def("add_node",[] (DirectedGraph  & ptr, DataType const & d, Size n_edges, Index n_end_index, NodeIndexandEdge const & pie) ->  Index {
		 return ptr.add_node(d, n_edges, n_end_index, pie); } )
		.def("add_node",[] (DirectedGraph  & ptr, DataType const & d, Size n_edges) ->  Index {
		 return ptr.add_node(d, n_edges); } )
		.def("has_parent",[] (DirectedGraph const & ptr, Index ni) ->  bool {
		 return ptr.has_parent(ni); } )
		.def("get_parent_index",[] (DirectedGraph const & ptr, Index ni) ->  Index {
		 return ptr.get_parent_index(ni); } )
		.def("get_parent_end_index",[] (DirectedGraph const & ptr, Index ni) ->  Index {
		 return ptr.get_parent_end_index(ni); } )
		.def("get_root_indexes",[] (DirectedGraph  & ptr) -> Indexes {
		 return ptr.get_root_indexes(); } )
		;

        py::class_<DirectedIterList, std::shared_ptr<DirectedIterList>>(m, "DirectedIterList")
		// ctors
		.def(py::init<>())
		// methods
		.def("transversal",[] (DirectedIterList  & ptr, AdjacencyListType & adj_list, Index start_n) {
		ptr.transversal(adj_list, start_n); } )
		.def("sub_graph_transversal",[] (DirectedIterList  & ptr, AdjacencyListType & adj_list, Index start_n, Index end_n) {
		ptr.sub_graph_transversal(adj_list, start_n, end_n); } )
		;
*/
    py::class_<DynamicEdges, std::shared_ptr<DynamicEdges>>(m, "DynamicEdges")
            ;

    py::class_<Edge, std::shared_ptr<Edge>>(m, "Edge")
            // ctors
            .def(py::init<Index,Index,Index,Index>())
                    // methods
            .def("partner",[] (Edge const & ptr, Index index) -> Index {
                return ptr.partner(index); } )
            .def("end_index",[] (Edge const & ptr, Index index) -> Index {
                return ptr.end_index(index); } )
            .def("to_str",[] (Edge const & ptr) -> String {
                return ptr.to_str(); } )
                    // operators
            .def(py::self == py::self)
                    // public attributes
            .def_readwrite("node_i", &Edge::node_i)
            .def_readwrite("node_j", &Edge::node_j)
            .def_readwrite("edge_i", &Edge::edge_i)
            .def_readwrite("edge_j", &Edge::edge_j)
            ;

    py::class_<FixedEdges, std::shared_ptr<FixedEdges>>(m, "FixedEdges")
            ;
/*
        py::class_<IterList, std::shared_ptr<IterList>>(m, "IterList")
		// ctors
		.def(py::init<>())
		// methods
		.def("begin",[] (IterList  & ptr) -> iterator {
		 return ptr.begin(); } )
		.def("end",[] (IterList  & ptr) -> iterator {
		 return ptr.end(); } )
		.def("begin",[] (IterList const & ptr) -> const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (IterList const & ptr) -> const_iterator {
		 return ptr.end(); } )
		.def("transversal",[] (IterList  & ptr, AdjacencyListType & adj_list, Index start_n) {
		ptr.transversal(adj_list, start_n); } )
		.def("path_transversal",[] (IterList  & ptr, AdjacencyListType & adj_list, Index start_n, Index end_n) {
		ptr.path_transversal(adj_list, start_n, end_n); } )
		;

        py::class_<Node, std::shared_ptr<Node>>(m, "Node")
		// ctors
		.def(py::init<DataType const &,Index const>())
		// methods
		.def("data",[] (Node const & ptr) ->  DataType const & {
		 return ptr.data(); } )
		.def("data",[] (Node  & ptr) ->  DataType & {
		 return ptr.data(); } )
		.def("index",[] (Node const & ptr) ->  Index {
		 return ptr.index(); } )
		;
*/
    py::class_<NodeIndexandEdge, std::shared_ptr<NodeIndexandEdge>>(m, "NodeIndexandEdge")
            // methods
            .def("to_str",[] (NodeIndexandEdge const & ptr) -> String {
                return ptr.to_str(); } )
                    // operators
            .def(py::self == py::self)
                    // public attributes
            .def_readwrite("node_index", &NodeIndexandEdge::node_index)
            .def_readwrite("edge_index", &NodeIndexandEdge::edge_index)
            ;
/*
        py::class_<NodeIndexandEdgeCompare, std::shared_ptr<NodeIndexandEdgeCompare>>(m, "NodeIndexandEdgeCompare")
		// operators
		.def(py::self () data_structure::NodeIndexandEdge const &)
		;

        py::class_<UndirectedGraph, std::shared_ptr<UndirectedGraph>>(m, "UndirectedGraph")
		// ctors
		.def(py::init<>())
		.def(py::init<UndirectedGraph const &>())
		;

        py::class_<VisitedNode, std::shared_ptr<VisitedNode>>(m, "VisitedNode")
		// methods
		.def("index_in_path",[] (VisitedNode  & ptr, Index i) -> bool {
		 return ptr.index_in_path(i); } )
		.def("path_length",[] (VisitedNode  & ptr) -> int {
		 return ptr.path_length(); } )
		// public attributes
		.def_readwrite("parent", &VisitedNode::parent)
		.def_readwrite("index", &VisitedNode::index)
		;

        py::class_<_Graph, std::shared_ptr<_Graph>>(m, "_Graph")
		// ctors
		.def(py::init<>())
		// methods
		.def("begin",[] (_Graph  & ptr) -> iterator {
		 return ptr.begin(); } )
		.def("end",[] (_Graph  & ptr) -> iterator {
		 return ptr.end(); } )
		.def("begin",[] (_Graph const & ptr) -> const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (_Graph const & ptr) -> const_iterator {
		 return ptr.end(); } )
		.def("setup_transversal",[] (_Graph  & ptr, Index start_n) {
		ptr.setup_transversal(start_n); } )
		.def("setup_path_transversal",[] (_Graph  & ptr, Index start_n, Index end_n) {
		ptr.setup_path_transversal(start_n, end_n); } )
		.def("add_node",[] (_Graph  & ptr, DataType const & d, Size n_edges) ->  Index {
		 return ptr.add_node(d, n_edges); } )
		.def("add_edge",[] (_Graph  & ptr, NodeIndexandEdge const & nie1, NodeIndexandEdge const & nie2) {
		ptr.add_edge(nie1, nie2); } )
		.def("remove_node",[] (_Graph  & ptr, Index ni) {
		ptr.remove_node(ni); } )
		.def("remove_edge",[] (_Graph  & ptr, NodeIndexandEdge const & nie1, NodeIndexandEdge const & nie2) {
		ptr.remove_edge(nie1, nie2); } )
		.def("get_num_nodes",[] (_Graph const & ptr) ->  size_t {
		 return ptr.get_num_nodes(); } )
		.def("get_num_edges",[] (_Graph const & ptr) ->  size_t {
		 return ptr.get_num_edges(); } )
		.def("get_node_edges",[] (_Graph const & ptr, Index ni) ->  std::vector<Edge const *> const & {
		 return ptr.get_node_edges(ni); } )
		.def("get_node",[] (_Graph const & ptr, Index ni) ->  Node<DataType> const & {
		 return ptr.get_node(ni); } )
		.def("get_node_data",[] (_Graph const & ptr, Index ni) ->  DataType const & {
		 return ptr.get_node_data(ni); } )
		.def("get_node_data",[] (_Graph  & ptr, Index ni) ->  DataType & {
		 return ptr.get_node_data(ni); } )
		.def("get_connected_node_info",[] (_Graph const & ptr, NodeIndexandEdge const & nei) ->  NodeIndexandEdge {
		 return ptr.get_connected_node_info(nei); } )
		.def("edge_between_nodes",[] (_Graph const & ptr, Index n1, Index n2) ->  bool {
		 return ptr.edge_between_nodes(n1, n2); } )
		.def("edge_index_empty",[] (_Graph const & ptr, Index ni, Index ei) ->  bool {
		 return ptr.edge_index_empty(ni, ei); } )
		;
*/

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// data_structure::graph
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // classes
/*
        py::class_<Graph, std::shared_ptr<Graph>>(m, "Graph")
		// ctors
		.def(py::init<>())
		// methods
		.def("begin",[] (Graph  & ptr) -> iterator {
		 return ptr.begin(); } )
		.def("end",[] (Graph  & ptr) -> iterator {
		 return ptr.end(); } )
		.def("begin",[] (Graph const & ptr) -> const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (Graph const & ptr) -> const_iterator {
		 return ptr.end(); } )
		.def("transverse",[] (Graph const & ptr, data_structure::graphGraphNodeOP<DataType> const &) -> iterator {
		 return ptr.transverse(&); } )
		.def("size",[] (Graph const & ptr) ->  size_t {
		 return ptr.size(); } )
		.def("get_node",[] (Graph const & ptr, int index) ->  GraphNodeOP<DataType> const & {
		 return ptr.get_node(index); } )
		.def("oldest_node",[] (Graph  & ptr) ->  GraphNodeOP<DataType> {
		 return ptr.oldest_node(); } )
		.def("increase_level",[] (Graph  & ptr) {
		ptr.increase_level(); } )
		.def("decrease_level",[] (Graph  & ptr) {
		ptr.decrease_level(); } )
		.def("nodes",[] (Graph const & ptr) ->  GraphNodeOPs<DataType> const & {
		 return ptr.nodes(); } )
		.def("connections",[] (Graph const & ptr) ->  GraphConnectionOPs<DataType> const & {
		 return ptr.connections(); } )
		.def("last_node",[] (Graph  & ptr) ->  GraphNodeOP<DataType> const & {
		 return ptr.last_node(); } )
		.def("level",[] (Graph  & ptr) ->  int {
		 return ptr.level(); } )
		.def("index",[] (Graph  & ptr) ->  int {
		 return ptr.index(); } )
		.def("index",[] (Graph  & ptr, int nindex) {
		ptr.index(nindex); } )
		;

        py::class_<GraphConnection, std::shared_ptr<GraphConnection>>(m, "GraphConnection")
		// ctors
		.def(py::init<GraphNodeOP<DataType>,GraphNodeOP<DataType>,int,int>())
		// methods
		.def("disconnect",[] (GraphConnection  & ptr) {
		ptr.disconnect(); } )
		.def("partner",[] (GraphConnection  & ptr, int i) -> GraphNodeOP<DataType> const & {
		 return ptr.partner(i); } )
		.def("end_index",[] (GraphConnection  & ptr, int n_index) -> int {
		 return ptr.end_index(n_index); } )
		.def("node_1",[] (GraphConnection  & ptr) ->  GraphNodeOP<DataType> const & {
		 return ptr.node_1(); } )
		.def("node_2",[] (GraphConnection  & ptr) ->  GraphNodeOP<DataType> const & {
		 return ptr.node_2(); } )
		.def("end_index_1",[] (GraphConnection  & ptr) ->  int {
		 return ptr.end_index_1(); } )
		.def("end_index_2",[] (GraphConnection  & ptr) ->  int {
		 return ptr.end_index_2(); } )
		;

        py::class_<GraphDynamic, std::shared_ptr<GraphDynamic>>(m, "GraphDynamic")
		// ctors
		.def(py::init<>())
		// methods
		.def("add_data",[] (GraphDynamic  & ptr, DataType const & data, int parent_index, int orphan) ->  int {
		 return ptr.add_data(data, parent_index, orphan); } )
		.def("connect",[] (GraphDynamic  & ptr, int i, int j) {
		ptr.connect(i, j); } )
		;

        py::class_<GraphIterator, std::shared_ptr<GraphIterator>>(m, "GraphIterator")
		// ctors
		.def(py::init<>())
		// operators
		.def(py::self ++ int)
		.def(py::self == py::self)
		.def(py::self != py::self)
		;

        py::class_<GraphNode, std::shared_ptr<GraphNode>>(m, "GraphNode")
		// ctors
		.def(py::init<int,int,size_t>())
		.def(py::init<DataType const &,int,int,size_t>())
		// methods
		.def("add_connection",[] (GraphNode  & ptr, data_structure::graphGraphConnectionOP<DataType> const &, int pos) {
		ptr.add_connection(&, pos); } )
		.def("remove_connection",[] (GraphNode  & ptr, data_structure::graphGraphConnectionOP<DataType> const & connection) {
		ptr.remove_connection(connection); } )
		.def("available_children_pos",[] (GraphNode const & ptr) ->  Ints {
		 return ptr.available_children_pos(); } )
		.def("available_pos",[] (GraphNode  & ptr, int pos) ->  int {
		 return ptr.available_pos(pos); } )
		.def("parent",[] (GraphNode  & ptr) ->  GraphNodeOP<DataType> {
		 return ptr.parent(); } )
		.def("parent_index",[] (GraphNode  & ptr) ->  int {
		 return ptr.parent_index(); } )
		.def("parent_end_index",[] (GraphNode  & ptr) ->  int {
		 return ptr.parent_end_index(); } )
		.def("unset_connections",[] (GraphNode  & ptr) {
		ptr.unset_connections(); } )
		.def("connected",[] (GraphNode  & ptr, data_structure::graphGraphNodeOP<DataType> const & n) ->  GraphConnectionOP<DataType> {
		 return ptr.connected(n); } )
		.def("index",[] (GraphNode const & ptr) ->  int {
		 return ptr.index(); } )
		.def("level",[] (GraphNode const & ptr) ->  int {
		 return ptr.level(); } )
		.def("data",[] (GraphNode const & ptr) ->  DataType const & {
		 return ptr.data(); } )
		.def("data",[] (GraphNode  & ptr) ->  DataType & {
		 return ptr.data(); } )
		.def("connections",[] (GraphNode const & ptr) ->  GraphConnectionOPs<DataType> const & {
		 return ptr.connections(); } )
		.def("index",[] (GraphNode  & ptr, int index) {
		ptr.index(index); } )
		;

        py::class_<GraphNodeCompare, std::shared_ptr<GraphNodeCompare>>(m, "GraphNodeCompare")
		// operators
		.def(py::self () GraphNodeOP<DataType> const &)
		;

        py::class_<GraphNodeDynamic, std::shared_ptr<GraphNodeDynamic>>(m, "GraphNodeDynamic")
		// ctors
		.def(py::init<DataType const &,int,int>())
		// methods
		.def("add_connection",[] (GraphNodeDynamic  & ptr, data_structure::graphGraphConnectionOP<DataType> const & connection, int pos) {
		ptr.add_connection(connection, pos); } )
		.def("remove_connection",[] (GraphNodeDynamic  & ptr, data_structure::graphGraphConnectionOP<DataType> const & connection) {
		ptr.remove_connection(connection); } )
		;

        py::class_<GraphNodeStatic, std::shared_ptr<GraphNodeStatic>>(m, "GraphNodeStatic")
		// ctors
		.def(py::init<DataType const &,int,int,int>())
		.def(py::init<GraphNode<DataType> const &>())
		// methods
		.def("add_connection",[] (GraphNodeStatic  & ptr, data_structure::graphGraphConnectionOP<DataType> const & connection, int pos) {
		ptr.add_connection(connection, pos); } )
		.def("remove_connection",[] (GraphNodeStatic  & ptr, data_structure::graphGraphConnectionOP<DataType> const & connection) {
		ptr.remove_connection(connection); } )
		;

        py::class_<GraphNodeType, std::shared_ptr<GraphNodeType>>(m, "GraphNodeType")
		;

        py::class_<GraphStatic, std::shared_ptr<GraphStatic>>(m, "GraphStatic")
		// ctors
		.def(py::init<>())
		.def(py::init<GraphStatic<DataType> const &>())
		// methods
		.def("add_data",[] (GraphStatic  & ptr, DataType const & data, int parent_index, int parent_pos, int child_pos, int n_children, int orphan, int index) ->  int {
		 return ptr.add_data(data, parent_index, parent_pos, child_pos, n_children, orphan, index); } )
		.def("connect",[] (GraphStatic  & ptr, int i, int j, int i_pos, int j_pos) {
		ptr.connect(i, j, i_pos, j_pos); } )
		.def("check_pos_is_valid",[] (GraphStatic  & ptr, data_structure::graphGraphNodeOP<DataType> const & n, int & pos) ->  int {
		 return ptr.check_pos_is_valid(n, pos); } )
		.def("get_available_pos",[] (GraphStatic  & ptr, data_structure::graphGraphNodeOP<DataType> const & n, int & pos) ->  Ints {
		 return ptr.get_available_pos(n, pos); } )
		.def("remove_node",[] (GraphStatic  & ptr, int pos) {
		ptr.remove_node(pos); } )
		.def("remove_level",[] (GraphStatic  & ptr, int level) {
		ptr.remove_level(level); } )
		;


*/

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// data_structure::tree
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // classes
/*
        py::class_<Tree, std::shared_ptr<Tree>>(m, "Tree")
		// ctors
		.def(py::init<>())
		// methods
		.def("begin",[] (Tree  & ptr) -> iterator {
		 return ptr.begin(); } )
		.def("end",[] (Tree  & ptr) -> iterator {
		 return ptr.end(); } )
		.def("begin",[] (Tree const & ptr) -> const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (Tree const & ptr) -> const_iterator {
		 return ptr.end(); } )
		.def("get_node",[] (Tree  & ptr, int index) ->  TreeNodeOP<DataType> const & {
		 return ptr.get_node(index); } )
		.def("remove_node",[] (Tree  & ptr, data_structure::treeTreeNodeOP<DataType> const & n) {
		ptr.remove_node(n); } )
		.def("remove_node",[] (Tree  & ptr, int pos) {
		ptr.remove_node(pos); } )
		.def("remove_level",[] (Tree  & ptr, int level) {
		ptr.remove_level(level); } )
		.def("increase_level",[] (Tree  & ptr) {
		ptr.increase_level(); } )
		.def("decrease_level",[] (Tree  & ptr) {
		ptr.decrease_level(); } )
		.def("size",[] (Tree  & ptr) ->  size_t {
		 return ptr.size(); } )
		.def("last_node",[] (Tree  & ptr) ->  TreeNodeOP<DataType> const & {
		 return ptr.last_node(); } )
		.def("level",[] (Tree  & ptr) ->  int {
		 return ptr.level(); } )
		;

        py::class_<TreeDynamic, std::shared_ptr<TreeDynamic>>(m, "TreeDynamic")
		// ctors
		.def(py::init<>())
		// methods
		.def("add_data",[] (TreeDynamic  & ptr, DataType const & data, int parent_index) ->  int {
		 return ptr.add_data(data, parent_index); } )
		;

        py::class_<TreeNode, std::shared_ptr<TreeNode>>(m, "TreeNode")
		// ctors
		.def(py::init<DataType const &,int,int,size_t>())
		// methods
		.def("add_child",[] (TreeNode  & ptr, data_structure::treeTreeNodeOP<DataType> const &, int pos) {
		ptr.add_child(&, pos); } )
		.def("remove_child",[] (TreeNode  & ptr, data_structure::treeTreeNodeOP<DataType> const &) {
		ptr.remove_child(&); } )
		.def("leaf",[] (TreeNode  & ptr) -> bool {
		 return ptr.leaf(); } )
		.def("available_children_pos",[] (TreeNode const & ptr) ->  Ints {
		 return ptr.available_children_pos(); } )
		.def("available_pos",[] (TreeNode  & ptr, int pos) ->  int {
		 return ptr.available_pos(pos); } )
		.def("unset_children",[] (TreeNode  & ptr) {
		ptr.unset_children(); } )
		.def("index",[] (TreeNode const & ptr) ->  int {
		 return ptr.index(); } )
		.def("level",[] (TreeNode const & ptr) ->  int {
		 return ptr.level(); } )
		.def("data",[] (TreeNode const & ptr) ->  DataType const & {
		 return ptr.data(); } )
		.def("data",[] (TreeNode  & ptr) ->  DataType & {
		 return ptr.data(); } )
		.def("children",[] (TreeNode  & ptr) ->  TreeNodeOPs<DataType> const & {
		 return ptr.children(); } )
		.def("parent",[] (TreeNode  & ptr) ->  TreeNodeOP<DataType> const & {
		 return ptr.parent(); } )
		.def("parent_index",[] (TreeNode  & ptr) ->  int {
		 return ptr.parent_index(); } )
		.def("parent_end_index",[] (TreeNode  & ptr) ->  int {
		 return ptr.parent_end_index(); } )
		.def("parent",[] (TreeNode  & ptr, data_structure::treeTreeNodeOP<DataType> const & p) {
		ptr.parent(p); } )
		.def("index",[] (TreeNode  & ptr, int index) {
		ptr.index(index); } )
		;

        py::class_<TreeNodeDynamic, std::shared_ptr<TreeNodeDynamic>>(m, "TreeNodeDynamic")
		// ctors
		.def(py::init<DataType const &,int,int>())
		// methods
		.def("add_child",[] (TreeNodeDynamic  & ptr, data_structure::treeTreeNodeOP<DataType> const & c, int pos) {
		ptr.add_child(c, pos); } )
		.def("remove_child",[] (TreeNodeDynamic  & ptr, data_structure::treeTreeNodeOP<DataType> const & child) {
		ptr.remove_child(child); } )
		.def("leaf",[] (TreeNodeDynamic  & ptr) ->  bool {
		 return ptr.leaf(); } )
		;

        py::class_<TreeNodeStatic, std::shared_ptr<TreeNodeStatic>>(m, "TreeNodeStatic")
		// ctors
		.def(py::init<DataType const &,int,int,int>())
		// methods
		.def("add_child",[] (TreeNodeStatic  & ptr, data_structure::treeTreeNodeOP<DataType> const & child, int pos) {
		ptr.add_child(child, pos); } )
		.def("remove_child",[] (TreeNodeStatic  & ptr, data_structure::treeTreeNodeOP<DataType> const & child) {
		ptr.remove_child(child); } )
		.def("leaf",[] (TreeNodeStatic  & ptr) ->  bool {
		 return ptr.leaf(); } )
		;

        py::class_<TreeStatic, std::shared_ptr<TreeStatic>>(m, "TreeStatic")
		// ctors
		.def(py::init<>())
		.def(py::init<TreeStatic<DataType> const &>())
		// methods
		.def("add_data",[] (TreeStatic  & ptr, DataType const & data, int n_children, int parent_index, int parent_child_index) ->  int {
		 return ptr.add_data(data, n_children, parent_index, parent_child_index); } )
		.def("get_available_pos",[] (TreeStatic  & ptr, data_structure::treeTreeNodeOP<DataType> const & n, int pos) ->  Ints {
		 return ptr.get_available_pos(n, pos); } )
		.def("index",[] (TreeStatic  & ptr, int n_index) {
		ptr.index(n_index); } )
		;

*/
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// util
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

    py::class_<secondary_structure::Basepair, std::shared_ptr<secondary_structure::Basepair>>(m, "SecondaryStructureBasepair")
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

    py::class_<secondary_structure::Chain, std::shared_ptr<secondary_structure::Chain>>(m, "SecondaryStructureChain")
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

    py::class_<secondary_structure::Motif, std::shared_ptr<secondary_structure::Motif>>(m, "SecondaryStructureMotif")
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

    py::class_<secondary_structure::Pose, std::shared_ptr<secondary_structure::Pose>>(m, "SecondaryStructurePose")
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

    py::class_<secondary_structure::RNAStructure, std::shared_ptr<secondary_structure::RNAStructure>>(m, "SecondaryStructureRNAStructure")
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

    py::class_<secondary_structure::Residue, std::shared_ptr<secondary_structure::Residue>>(m, "SecondaryStructureResidue")
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
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// eternabot
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // classes

    py::class_<ABasicTest, std::shared_ptr<ABasicTest>>(m, "ABasicTest")
            // ctors
            .def(py::init<>())
                    // methods
            .def("score",[] (ABasicTest  & ptr, FeaturesOP const & features) ->  float {
                return ptr.score(features); } )
                    // inherited methods

            .def("mean",[] (Strategy const & ptr) ->  float {
                return ptr.mean(); } )
            .def("stdev",[] (Strategy const & ptr) ->  float {
                return ptr.stdev(); } )
            ;

    py::class_<BerexTest, std::shared_ptr<BerexTest>>(m, "BerexTest")
            // ctors
            .def(py::init<>())
                    // methods
            .def("score",[] (BerexTest  & ptr, FeaturesOP const & features) -> float {
                return ptr.score(features); } )
                    // inherited methods
            .def("mean",[] (Strategy const & ptr) ->  float {
                return ptr.mean(); } )
            .def("stdev",[] (Strategy const & ptr) ->  float {
                return ptr.stdev(); } )
            ;

    py::class_<CleanPlotStackCapsandSafeGC, std::shared_ptr<CleanPlotStackCapsandSafeGC>>(m, "CleanPlotStackCapsandSafeGC")
            // ctors
            .def(py::init<>())
                    // methods
            .def("score",[] (CleanPlotStackCapsandSafeGC  & ptr, FeaturesOP const & features) -> float {
                return ptr.score(features); } )
                    // inherited methods
            .def("mean",[] (Strategy const & ptr) ->  float {
                return ptr.mean(); } )
            .def("stdev",[] (Strategy const & ptr) ->  float {
                return ptr.stdev(); } )
            ;

    py::class_<DirectionofGCPairsinMultiLoops, std::shared_ptr<DirectionofGCPairsinMultiLoops>>(m, "DirectionofGCPairsinMultiLoops")
            // ctors
            .def(py::init<>())
                    // methods
            .def("score",[] (DirectionofGCPairsinMultiLoops  & ptr, FeaturesOP const & features) -> float {
                return ptr.score(features); } )
                    // inherited methods
            .def("mean",[] (Strategy const & ptr) ->  float {
                return ptr.mean(); } )
            .def("stdev",[] (Strategy const & ptr) ->  float {
                return ptr.stdev(); } )
            ;

    py::class_<FeatureGenerator, std::shared_ptr<FeatureGenerator>>(m, "FeatureGenerator")
            // ctors
            .def(py::init<>())
                    // methods
            .def("get_features",[] (FeatureGenerator  & ptr, secondary_structure::PoseOP const & p) -> FeaturesOP {
                return ptr.get_features(p); } )
            .def("update_features",[] (FeatureGenerator  & ptr, FeaturesOP & features, secondary_structure::PoseOP const & p) {
                ptr.update_features(features, p); } )
            ;

    py::class_<Features, std::shared_ptr<Features>>(m, "Features")
            // ctors
            .def(py::init<>())
                    // public attributes
            .def_readwrite("length", &Features::length)
            .def_readwrite("a_count", &Features::a_count)
            .def_readwrite("c_count", &Features::c_count)
            .def_readwrite("g_count", &Features::g_count)
            .def_readwrite("u_count", &Features::u_count)
            .def_readwrite("gu", &Features::gu)
            .def_readwrite("gc", &Features::gc)
            .def_readwrite("ua", &Features::ua)
            .def_readwrite("meltpoint", &Features::meltpoint)
            .def_readwrite("fe", &Features::fe)
            .def_readwrite("dotplot", &Features::dotplot)
            .def_readwrite("pairmap", &Features::pairmap)
            .def_readwrite("helices", &Features::helices)
            ;

    py::class_<NumofYellowNucleotidesperLengthofString, std::shared_ptr<NumofYellowNucleotidesperLengthofString>>(m, "NumofYellowNucleotidesperLengthofString")
            // ctors
            .def(py::init<>())
                    // methods
            .def("score",[] (NumofYellowNucleotidesperLengthofString  & ptr, FeaturesOP const & features) -> float {
                return ptr.score(features); } )
                    // inherited methods
            .def("mean",[] (Strategy const & ptr) ->  float {
                return ptr.mean(); } )
            .def("stdev",[] (Strategy const & ptr) ->  float {
                return ptr.stdev(); } )
            ;
/*
        py::class_<Parameters, std::shared_ptr<Parameters>>(m, "Parameters")
		;
*/
    py::class_<eternabot::Scorer, std::shared_ptr<eternabot::Scorer>>(m, "EternabotScorer")
            // ctors
            .def(py::init<>())
                    // methods
            .def("setup",[] (eternabot::Scorer  & ptr, secondary_structure::PoseOP const & p) {
                ptr.setup(p); } )
            .def("score_secondary_structure",[] (eternabot::Scorer  & ptr, secondary_structure::PoseOP const & p) -> float {
                return ptr.score_secondary_structure(p); } )
            .def("scores",[] (eternabot::Scorer  & ptr) -> Floats const & {
                return ptr.scores(); } )
            ;

    py::class_<SequenceDesigner, std::shared_ptr<SequenceDesigner>>(m, "SequenceDesigner")
            // ctors
            .def(py::init<>())
                    // methods
            .def("setup",[] (SequenceDesigner  & ptr) {
                ptr.setup(); } )
            .def("design",[] (SequenceDesigner  & ptr, secondary_structure::PoseOP const & p) -> SequenceDesignerResultOPs const &  {
                return ptr.design(p); } )
            .def("get_int_option",[] (SequenceDesigner  & ptr, String const & name) ->  float {
                return ptr.get_int_option(name); } )
            .def("get_float_option",[] (SequenceDesigner  & ptr, String const & name) ->  float {
                return ptr.get_float_option(name); } )
            .def("get_string_option",[] (SequenceDesigner  & ptr, String const & name) ->  String {
                return ptr.get_string_option(name); } )
            .def("get_bool_option",[] (SequenceDesigner  & ptr, String const & name) ->  bool {
                return ptr.get_bool_option(name); } )
            .def("has_option",[] (SequenceDesigner  & ptr, String const & name) ->  bool {
                return ptr.has_option(name); } )
            .def("set_option_value",[] (SequenceDesigner  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (SequenceDesigner  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (SequenceDesigner  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (SequenceDesigner  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )
            ;

    py::class_<SequenceDesignerResult, std::shared_ptr<SequenceDesignerResult>>(m, "SequenceDesignerResult")
            // ctors
            .def(py::init<String const &,float>())
                    // public attributes
            .def_readwrite("sequence", &SequenceDesignerResult::sequence)
            .def_readwrite("score", &SequenceDesignerResult::score)
            ;
/*
        py::class_<Strategy, std::shared_ptr<Strategy>>(m, "Strategy")
		// ctors
		.def(py::init<>())
		// methods
		.def("score",[] (Strategy  & ptr, FeaturesOP const &) -> float {
		 return ptr.score(&); } )
		.def("mean",[] (Strategy const & ptr) ->  float {
		 return ptr.mean(); } )
		.def("stdev",[] (Strategy const & ptr) ->  float {
		 return ptr.stdev(); } )
		;
        py::class_<sequence_designer_result_less_than_key, std::shared_ptr<sequence_designer_result_less_than_key>>(m, "sequence_designer_result_less_than_key")
		// operators
		.def(py::self () eternabot::SequenceDesignerResultOP)
		;

*/
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// structure
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

    m.def("are_chains_equal", [] (structure::ChainOP const & c1, structure::ChainOP const & c2, int check_uuids) -> bool {
              return structure::are_chains_equal(c1, c2, check_uuids); },
          py::arg("c1"),
          py::arg("c2"),
          py::arg("check_uuids") = 1
    );

    m.def("are_residues_equal", [] (structure::ResidueOP const & r1, structure::ResidueOP const & r2, int check_uuids) -> bool {
              return structure::are_residues_equal(r1, r2, check_uuids); },
          py::arg("r1"),
          py::arg("r2"),
          py::arg("check_uuids") = 1
    );

    m.def("are_structures_equal", [] (structure::StructureOP const & s1, structure::StructureOP const & s2, int check_uuids) -> bool {
              return structure::are_structures_equal(s1, s2, check_uuids); },
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

    m.def("close_chain", [] (structure::ChainOP chain) {
              structure::close_chain(chain); },
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

    m.def("connect_residues_into_chains", [] (structure::ResidueOPs & residues, structure::ChainOPs & chains) {
              structure::connect_residues_into_chains(residues, chains); },
          py::arg("residues"),
          py::arg("chains")
    );

    m.def("create_coord_system", [] (AtomOPs const & atoms) -> math::Matrix {
              return create_coord_system(atoms); },
          py::arg("atoms")
    );

    m.def("end_from_basepairs", [] (structure::StructureOP const & s, structure::BasepairOPs const & bps) {//-> std::shared_ptr<BasepairOPs> {
              return *structure::end_from_basepairs(s, bps); },
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

    m.def("get_res_ref_frame", [] (structure::ResidueOP r) -> math::Matrix {
              return structure::get_res_ref_frame(r); },
          py::arg("r")
    );

    m.def("gu_bp", [] (structure::BasepairOP const & bp) -> bool {
              return structure::gu_bp(bp); },
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

    m.def("replace_missing_phosphate_backbone", [] (structure::ResidueOP r, structure::ResidueOP r_template) {
              structure::replace_missing_phosphate_backbone(r, r_template); },
          py::arg("r"),
          py::arg("r_template")
    );

    m.def("str_to_basepairstate", [] (String const & s) -> BasepairState {
              return str_to_basepairstate(s); },
          py::arg("s")
    );

    m.def("subselect_basepairs_with_res", [] (structure::ResidueOPs const & res, structure::BasepairOPs const & all_bps) {
              return *structure::subselect_basepairs_with_res(res, all_bps); },
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

    m.def("wc_bp", [] (structure::BasepairOP const & bp) -> bool {
              return structure::wc_bp(bp); },
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

    py::class_<structure::Basepair, std::shared_ptr<structure::Basepair>>(m, "StructureBasepair")
            // ctors
            .def(py::init<>())
            .def(py::init<structure::ResidueOP &,structure::ResidueOP &,math::Matrix const &,String const &>())
                    // methods
            .def("copy",[] (structure::Basepair  & ptr) -> structure::Basepair {
                return ptr.copy(); } )
            .def("move",[] (structure::Basepair  & ptr, math::Point const & p) {
                ptr.move(p); } )
            .def("transform",[] (structure::Basepair  & ptr, math::Matrix const & r, math::Vector const & t, math::Point & dummy) {
                ptr.transform(r, t, dummy); } )
            .def("transform",[] (structure::Basepair  & ptr, math::Matrix const & r, math::Vector const & t) {
                ptr.transform(r, t); } )
            .def("state",[] (structure::Basepair  & ptr) ->  structure::BasepairStateOP const & {
                return ptr.state(); } )
            .def("diff",[] (structure::Basepair  & ptr, structure::Basepair & bp) ->  float {
                return ptr.diff(bp); } )
            .def("residues",[] (structure::Basepair const & ptr) ->  structure::ResidueOPs const {
                return ptr.residues(); } )
            .def("partner",[] (structure::Basepair  & ptr, structure::ResidueOP const & res) ->  structure::ResidueOP const & {
                return ptr.partner(res); } )
            .def("name",[] (structure::Basepair const & ptr) ->  String const {
                return ptr.name(); } )
            .def("flip",[] (structure::Basepair  & ptr) {
                ptr.flip(); } )
            .def("r",[] (structure::Basepair const & ptr) ->  math::Matrix const & {
                return ptr.r(); } )
            .def("d",[] (structure::Basepair const & ptr) ->  math::Point const {
                return ptr.d(); } )
            .def("uuid",[] (structure::Basepair const & ptr) ->  util::Uuid const & {
                return ptr.uuid(); } )
            .def("res1",[] (structure::Basepair const & ptr) ->  structure::ResidueOP {
                return ptr.res1(); } )
            .def("res2",[] (structure::Basepair const & ptr) ->  structure::ResidueOP {
                return ptr.res2(); } )
            .def("bp_type",[] (structure::Basepair const & ptr) ->  String const & {
                return ptr.bp_type(); } )
            .def("flipped",[] (structure::Basepair const & ptr) ->  int const {
                return ptr.flipped(); } )
            .def("atoms",[] (structure::Basepair const & ptr) ->  AtomOPs const & {
                return ptr.atoms(); } )
            .def("r",[] (structure::Basepair  & ptr, math::Matrix const & nr) {
                ptr.r(nr); } )
            .def("uuid",[] (structure::Basepair  & ptr, util::Uuid const & nuuid) {
                ptr.uuid(nuuid); } )
            .def("flipped",[] (structure::Basepair  & ptr, int const & nflipped) {
                ptr.flipped(nflipped); } )
            .def("res1",[] (structure::Basepair  & ptr, structure::ResidueOP const & nres1) {
                ptr.res1(nres1); } )
            .def("res2",[] (structure::Basepair  & ptr, structure::ResidueOP const & nres2) {
                ptr.res2(nres2); } )
            .def("bp_type",[] (structure::Basepair  & ptr, String const & nbp_type) {
                ptr.bp_type(nbp_type); } )
            .def("to_str",[] (structure::Basepair const & ptr) -> String const {
                return ptr.to_str(); } )
            .def("to_pdb_str",[] (structure::Basepair const & ptr) -> String const {
                return ptr.to_pdb_str(); } )
            .def("to_pdb",[] (structure::Basepair const & ptr, String const fname) {
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

    py::enum_<BeadType>(m, "BeadType")
        .value("PHOS", BeadType::PHOS)
        .value("SUGAR", BeadType::SUGAR)
        .value("BASE", BeadType::BASE)
            ;

    py::class_<CIFParser, std::shared_ptr<CIFParser>>(m, "CIFParser")
            // ctors
            .def(py::init<>())
                    // methods
            .def("parse",[] (CIFParser  & ptr, String const & pdb_file, int protein, int rna, int others) -> structure::ResidueOPs const & {
                return ptr.parse(pdb_file, protein, rna, others); } )
            ;

    py::class_<structure::Chain, std::shared_ptr<structure::Chain>>(m, "StructureChain")
            // ctors
            .def(py::init<>())
            .def(py::init<structure::ResidueOPs const &>())
            .def(py::init<structure::Chain const &>())
            .def(py::init<String const &,ResidueTypeSet const &>())
                    // methods
            .def("length",[] (structure::Chain const & ptr) ->  int {
                return ptr.length(); } )
            .def("first",[] (structure::Chain  & ptr) ->  structure::ResidueOP const & {
                return ptr.first(); } )
            .def("last",[] (structure::Chain  & ptr) ->  structure::ResidueOP const & {
                return ptr.last(); } )
            .def("subchain",[] (structure::Chain  & ptr, int start, int end) ->  structure::ChainOP {
                return ptr.subchain(start, end); } )
            .def("subchain",[] (structure::Chain  & ptr, structure::ResidueOP const & r1, structure::ResidueOP const & r2) ->  structure::ChainOP {
                return ptr.subchain(r1, r2); } )
            .def("to_str",[] (structure::Chain const & ptr) -> String {
                return ptr.to_str(); } )
            .def("to_pdb_str",[] (structure::Chain const & ptr, int & acount, int rnum, String const & chain_id) -> String {
                return ptr.to_pdb_str(acount, rnum, chain_id); } )
            .def("to_pdb_str",[] (structure::Chain const & ptr, int & acount) ->  String {
                return ptr.to_pdb_str(acount); } )
            .def("to_pdb",[] (structure::Chain const & ptr, String const fname, int rnum, String const & chain_id) {
                ptr.to_pdb(fname, rnum, chain_id); } )
            .def("to_pdb",[] (structure::Chain  & ptr, String const & fname) {
                ptr.to_pdb(fname); } )
            .def("residues",[] (structure::Chain  & ptr) ->  structure::ResidueOPs & {
                return ptr.residues(); } )
            ;

    py::class_<PDBParser, std::shared_ptr<PDBParser>>(m, "PDBParser")
            // ctors
            .def(py::init<>())
                    // methods
            .def("parse",[] (PDBParser  & ptr, String const & pdb_file, int protein, int rna, int others) -> structure::ResidueOPs const & {
                return ptr.parse(pdb_file, protein, rna, others); } )
            ;

    py::class_<structure::RNAStructure, std::shared_ptr<structure::RNAStructure>>(m, "StructureRNAStructure")
            // ctors
            .def(py::init<>())
            .def(py::init<structure::StructureOP const &,structure::BasepairOPs const &,structure::BasepairOPs const &>())
            .def(py::init<structure::RNAStructure const &>())
                    // methods
            .def("get_basepair",[] (structure::RNAStructure  & ptr, String const & str) -> structure::BasepairOPs {
                return ptr.get_basepair(str); } )
            .def("get_basepair",[] (structure::RNAStructure  & ptr, util::Uuid const & id) -> structure::BasepairOPs {
                return ptr.get_basepair(id); } )
            .def("get_basepair",[] (structure::RNAStructure  & ptr, structure::ResidueOP const & r1, structure::ResidueOP const & r2) -> structure::BasepairOPs {
                return ptr.get_basepair(r1, r2); } )
            .def("get_basepair",[] (structure::RNAStructure  & ptr, util::Uuid const & id1, util::Uuid const & id2) -> structure::BasepairOPs {
                return ptr.get_basepair(id1, id2); } )
            .def("get_beads",[] (structure::RNAStructure  & ptr, structure::BasepairOPs const & bps) -> Beads const & {
                return ptr.get_beads(bps); } )
            .def("get_beads",[] (structure::RNAStructure  & ptr, structure::BasepairOP const & bps) -> Beads const & {
                return ptr.get_beads(bps); } )
            .def("get_beads",[] (structure::RNAStructure  & ptr) ->  Beads const & {
                return ptr.get_beads(); } )
            .def("get_end_index",[] (structure::RNAStructure  & ptr, structure::BasepairOP const & bp) -> int {
                return ptr.get_end_index(bp); } )
            .def("get_end_index",[] (structure::RNAStructure  & ptr, String const & str) -> int {
                return ptr.get_end_index(str); } )
            .def("to_pdb_str",[] (structure::RNAStructure  & ptr, int rnumber, int close_chains) -> String const {
                return ptr.to_pdb_str(rnumber, close_chains); } )
            .def("to_pdb",[] (structure::RNAStructure  & ptr, String const fname, int renumber, int close_chains, int conect_statements) {
                ptr.to_pdb(fname, renumber, close_chains, conect_statements); } )
            .def("atoms",[] (structure::RNAStructure const & ptr) ->  AtomOPs const {
                return ptr.atoms(); } )
            .def("residues",[] (structure::RNAStructure const & ptr) ->  structure::ResidueOPs const {
                return ptr.residues(); } )
            .def("chains",[] (structure::RNAStructure  & ptr) -> structure::ChainOPs const & {
                return ptr.chains(); } )
            .def("get_residue",[] (structure::RNAStructure  & ptr, int num, String const & chain_id, String const & i_code) ->  structure::ResidueOP const {
                return ptr.get_residue(num, chain_id, i_code); } )
            .def("get_residue",[] (structure::RNAStructure  & ptr, util::Uuid const & uuid) ->  structure::ResidueOP const {
                return ptr.get_residue(uuid); } )
            .def("ends",[] (structure::RNAStructure const & ptr) ->  structure::BasepairOPs const & {
                return ptr.ends(); } )
            .def("end_ids",[] (structure::RNAStructure const & ptr) ->  Strings const & {
                return ptr.end_ids(); } )
            .def("name",[] (structure::RNAStructure const & ptr) ->  String const & {
                return ptr.name(); } )
            .def("path",[] (structure::RNAStructure const & ptr) ->  String const & {
                return ptr.path(); } )
            .def("basepairs",[] (structure::RNAStructure const & ptr) ->  structure::BasepairOPs const & {
                return ptr.basepairs(); } )
            .def("UNSAFE_basepairs",[] (structure::RNAStructure  & ptr) ->  structure::BasepairOPs & {
                return ptr.UNSAFE_basepairs(); } )
            .def("beads",[] (structure::RNAStructure const & ptr) ->  Beads const & {
                return ptr.beads(); } )
            .def("score",[] (structure::RNAStructure const & ptr) ->  float const & {
                return ptr.score(); } )
            .def("protein_beads",[] (structure::RNAStructure  & ptr) ->  Beads const & {
                return ptr.protein_beads(); } )
            .def("name",[] (structure::RNAStructure  & ptr, String const & nname) {
                ptr.name(nname); } )
            .def("path",[] (structure::RNAStructure  & ptr, String const & npath) {
                ptr.path(npath); } )
            .def("score",[] (structure::RNAStructure  & ptr, float const & nscore) {
                ptr.score(nscore); } )
            .def("end_ids",[] (structure::RNAStructure  & ptr, Strings const & end_ids) {
                ptr.end_ids(end_ids); } )
            .def("ends",[] (structure::RNAStructure  & ptr, structure::BasepairOPs const & ends) {
                ptr.ends(ends); } )
            .def("protein_beads",[] (structure::RNAStructure  & ptr, Beads const & beads) {
                ptr.protein_beads(beads); } )
            ;

    py::class_<structure::Residue, std::shared_ptr<structure::Residue>>(m, "StructureResidue")
            // ctors
            .def(py::init<structure::ResidueType const &,String const &,int const &,String const &,String const &>())
            .def(py::init<structure::Residue const &>())
            .def(py::init<String const &,structure::ResidueTypeSet const &>())
                    // methods
            .def("setup_atoms",[] (structure::Residue  & ptr, AtomOPs & atoms) {
                ptr.setup_atoms(atoms); } )
            .def("get_atom",[] (structure::Residue const & ptr, String const & name) ->  AtomOP const & {
                return ptr.get_atom(name); } )
            .def("connected_to",[] (structure::Residue const & ptr, structure::Residue const & res, float cutoff) ->  int {
                return ptr.connected_to(res, cutoff); } )
            .def("new_uuid",[] (structure::Residue  & ptr) {
                ptr.new_uuid(); } )
            .def("get_beads",[] (structure::Residue const & ptr) -> Beads {
                return ptr.get_beads(); } )
            .def("to_str",[] (structure::Residue const & ptr) -> String {
                return ptr.to_str(); } )
            .def("to_pdb_str",[] (structure::Residue  & ptr, int & acount) ->  String {
                return ptr.to_pdb_str(acount); } )
            .def("to_pdb_str",[] (structure::Residue const & ptr, int & acount, int rnum, String const & chain_id) -> String {
                return ptr.to_pdb_str(acount, rnum, chain_id); } )
            .def("to_pdb",[] (structure::Residue  & ptr, String const fname) {
                ptr.to_pdb(fname); } )
            .def("move",[] (structure::Residue  & ptr, math::Point const & p) {
                ptr.move(p); } )
            .def("transform",[] (structure::Residue  & ptr, math::Transform const & t) {
                ptr.transform(t); } )
            .def("num",[] (structure::Residue  & ptr, int nnum) {
                ptr.num(nnum); } )
            .def("chain_id",[] (structure::Residue  & ptr, String const & nchain_id) {
                ptr.chain_id(nchain_id); } )
            .def("uuid",[] (structure::Residue  & ptr, util::Uuid const & uuid) {
                ptr.uuid(uuid); } )
            .def("center",[] (structure::Residue const & ptr) ->  math::Point {
                return ptr.center(); } )
            .def("name",[] (structure::Residue const & ptr) ->  String const & {
                return ptr.name(); } )
            .def("chain_id",[] (structure::Residue const & ptr) ->  String const & {
                return ptr.chain_id(); } )
            .def("i_code",[] (structure::Residue const & ptr) ->  String const & {
                return ptr.i_code(); } )
            .def("num",[] (structure::Residue const & ptr) ->  int {
                return ptr.num(); } )
            .def("short_name",[] (structure::Residue const & ptr) ->  String {
                return ptr.short_name(); } )
            .def("atoms",[] (structure::Residue const & ptr) ->  AtomOPs const & {
                return ptr.atoms(); } )
            .def("uuid",[] (structure::Residue const & ptr) ->  util::Uuid const & {
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

    py::enum_<SetType>(m, "SetType")
            .value("RNA", SetType::RNA)
            .value("PROTEIN", SetType::PROTEIN)
            .value("UNKNOWN", SetType::UNKNOWN)
            ;

    py::class_<structure::Structure, std::shared_ptr<structure::Structure>>(m, "StructureStructure")
            // ctors
            .def(py::init<>())
            .def(py::init<structure::ChainOPs const &>())
            .def(py::init<String const &>())
            .def(py::init<structure::Structure const &>())
            .def(py::init<String const &,ResidueTypeSet const &>())
                    // methods
            .def("renumber",[] (structure::Structure  & ptr) {
                ptr.renumber(); } )
            .def("get_beads",[] (structure::Structure  & ptr, structure::ResidueOPs const & excluded_res) ->  Beads {
                return ptr.get_beads(excluded_res); } )
            .def("get_beads",[] (structure::Structure  & ptr) ->  Beads {
                return ptr.get_beads(); } )
            .def("get_residue",[] (structure::Structure  & ptr, int const & num, String const & chain_id, String const & i_code) -> structure::ResidueOP {
                return ptr.get_residue(num, chain_id, i_code); } )
            .def("get_residue",[] (structure::Structure  & ptr, int const & num, String const & chain_id, String const & i_code) -> structure::ResidueOP {
                return ptr.get_residue(num, chain_id, i_code); } )
            .def("residues",[] (structure::Structure const & ptr) -> structure::ResidueOPs const {
                return ptr.residues(); } )
            .def("atoms",[] (structure::Structure  & ptr) ->  AtomOPs const {
                return ptr.atoms(); } )
            .def("move",[] (structure::Structure  & ptr, math::Point const & p) {
                ptr.move(p); } )
            .def("transform",[] (structure::Structure  & ptr, math::Transform const & t) {
                ptr.transform(t); } )
            .def("to_pdb_str",[] (structure::Structure  & ptr, int renumber, int conect_statements) -> String {
                return ptr.to_pdb_str(renumber, conect_statements); } )
            .def("to_str",[] (structure::Structure  & ptr) -> String {
                return ptr.to_str(); } )
            .def("to_pdb",[] (structure::Structure  & ptr, String const fname, int renumber, int conect_statements) {
                ptr.to_pdb(fname, renumber, conect_statements); } )
            .def("chains",[] (structure::Structure const & ptr) ->  structure::ChainOPs const & {
                return ptr.chains(); } )
            ;
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

    py::class_<motif::Motif, std::shared_ptr<motif::Motif>>(m, "MotifMotif")
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

    py::class_<motif::Pose, std::shared_ptr<motif::Pose>>(m, "MotifPose")
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
            .def("get_beads",[] (motif::Pose  & ptr, structure::BasepairOPs const & bps) -> Beads const & {
                return ptr.get_beads(bps); } )
            .def("get_beads",[] (motif::Pose  & ptr, structure::BasepairOP const & bps) -> Beads const & {
                return ptr.get_beads(bps); } )
            .def("get_beads",[] (motif::Pose  & ptr) ->  Beads const & {
                return ptr.get_beads(); } )
            .def("get_end_index",[] (motif::Pose  & ptr, structure::BasepairOP const & bp) -> int {
                return ptr.get_end_index(bp); } )
            .def("get_end_index",[] (motif::Pose  & ptr, String const & str) -> int {
                return ptr.get_end_index(str); } )
            .def("to_pdb_str",[] (motif::Pose  & ptr, int rnumber, int close_chains) -> String const {
                return ptr.to_pdb_str(rnumber, close_chains); } )
            .def("to_pdb",[] (motif::Pose  & ptr, String const fname, int renumber, int close_chains, int conect_statements) {
                ptr.to_pdb(fname, renumber, close_chains, conect_statements); } )
            .def("atoms",[] (motif::Pose const & ptr) ->  AtomOPs const {
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
            .def("beads",[] (motif::Pose const & ptr) ->  Beads const & {
                return ptr.beads(); } )
            .def("score",[] (motif::Pose const & ptr) ->  float const & {
                return ptr.score(); } )
            .def("protein_beads",[] (motif::Pose  & ptr) ->  Beads const & {
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
            .def("protein_beads",[] (motif::Pose  & ptr, Beads const & beads) {
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

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// motif_tools
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // classes

    py::class_<Pair, std::shared_ptr<Pair>>(m, "Pair")
            // ctors
            .def(py::init<structure::ResidueOP const &,structure::ResidueOP const &,int>())
                    // methods
            .def("contains",[] (Pair const & ptr, structure::ResidueOP const & res) -> int  {
                return ptr.contains(res); } )
                    // public attributes
            .def_readwrite("res1", &Pair::res1)
            .def_readwrite("res2", &Pair::res2)
            .def_readwrite("dist", &Pair::dist)
            ;

    py::class_<PairSearch, std::shared_ptr<PairSearch>>(m, "PairSearch")
            // ctors
            .def(py::init<>())
                    // methods
            .def("search",[] (PairSearch  & ptr, structure::ResidueOPs const & res, PairOPs pairs, PairOPs end_pairs) -> PairSearchNodes const & {
                return ptr.search(res, pairs, end_pairs); } )
            ;

    py::class_<PairSearchNode, std::shared_ptr<PairSearchNode>>(m, "PairSearchNode")
            // ctors
            .def(py::init<motif_tools::PairOPs>())
                    // methods
            .def("contains",[] (PairSearchNode const & ptr, structure::ResidueOP const & res) ->  int {
                return ptr.contains(res); } )
            .def("contains_pair",[] (PairSearchNode  & ptr, PairOP pair) ->  int {
                return ptr.contains_pair(pair); } )
                    // public attributes
            .def_readwrite("pairs", &PairSearchNode::pairs)
            .def_readwrite("score", &PairSearchNode::score)
            ;
/*
        py::class_<PairSearchNodeCompare, std::shared_ptr<PairSearchNodeCompare>>(m, "PairSearchNodeCompare")
		// operators
		.def(py::self () PairSearchNode const &)
		;
*/
    py::class_<Segmenter, std::shared_ptr<Segmenter>>(m, "Segmenter")
            // ctors
            .def(py::init<>())
                    // methods
            .def("apply",[] (Segmenter  & ptr, structure::RNAStructureOP const & m,structure::BasepairOPs const & bps) -> SegmentsOP {
                return ptr.apply(m, bps); } )
            ;

    py::class_<Segments, std::shared_ptr<Segments>>(m, "Segments")
            // ctors
            .def(py::init<motif::MotifOP const &,motif::MotifOP const &>())
                    // public attributes
            .def_readwrite("removed", &Segments::removed)
            .def_readwrite("remaining", &Segments::remaining)
            ;


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// resources
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // free functions

    m.def("build_sqlite_library", [] (String const & path, std::vector<Strings>const & data, Strings const & keys, String const & primary_key) {
              build_sqlite_library(path, data, keys, primary_key); },
          py::arg("path"),
          py::arg("data"),
          py::arg("keys"),
          py::arg("primary_key")
    );

    m.def("get_motif_from_resource_manager", [] (String const & name, String const & end_id, String const & end_name) ->  motif::MotifOP {
              return get_motif_from_resource_manager(name, end_id, end_name); },
          py::arg("name") = dummy_name,
          py::arg("end_id") = dummy_end_id,
          py::arg("end_name") = dummy_name
    );

    m.def("sqlite3_escape", [] (String & unescaped_string) {
              sqlite3_escape(unescaped_string); },
          py::arg("unescaped_string")
    );
    // classes

    py::class_<AddedMotifLibrary, std::shared_ptr<AddedMotifLibrary>>(m, "AddedMotifLibrary")
            // ctors
            .def(py::init<>())
                    // methods
            .def("add_motif",[] (AddedMotifLibrary  & ptr, motif::MotifOP const & m) {
                ptr.add_motif(m); } )
            .def("get",[] (AddedMotifLibrary  & ptr, String const & name, String const & end_id, String const & end_name) -> motif::MotifOP {
                return ptr.get(name, end_id, end_name); } )
            .def("get_multi",[] (AddedMotifLibrary  & ptr, String const & name, String const & end_id, String const & end_name) -> motif::MotifOPs {
                return ptr.get_multi(name, end_id, end_name); } )
            .def("contains",[] (AddedMotifLibrary  & ptr, String const & name, String const & end_id, String const & end_name) -> int {
                return ptr.contains(name, end_id, end_name); } )
            ;

    py::class_<Manager, std::unique_ptr<Manager,py::nodelete>>(m, "Manager")
            .def(py::init([](){
                return std::unique_ptr<Manager,py::nodelete>(&Manager::instance());
            }))
                    //.def(py::init([](){return Manager::instance();}))
                    //.def("__init__",[]() {return Manager::instance();})
                    //.def(py::init([](){return Manager::instance();}))
                    // methods
                    //.def(py::init(),py::return_value_policy::automatic_reference)
            .def("instance",[] (Manager  & ptr) -> Manager & {
                return ptr.instance(); } )
            .def("bp_step",[] (Manager  & ptr, String const & end_id) -> motif::MotifOP {
                return ptr.bp_step(end_id); } )
            .def("bp_step_state",[] (Manager  & ptr, String const & end_id) -> motif::MotifStateOP {
                return ptr.bp_step_state(end_id); } )
            .def("motif",[] (Manager  & ptr, String const & name, String const & end_id, String const & end_name) -> motif::MotifOP {
                return ptr.motif(name, end_id, end_name); } )
            .def("motif_state",[] (Manager  & ptr, String const & name, String const & end_id, String const & end_name) -> motif::MotifStateOP {
                return ptr.motif_state(name, end_id, end_name); } )
            .def("motif_state_ensemble",[] (Manager  & ptr, String const & name) -> motif::MotifStateEnsembleOP {
                return ptr.motif_state_ensemble(name); } )
            .def("get_structure",[] (Manager  & ptr, String const & path, String name, int force_num_chains) -> structure::RNAStructureOP {
                return ptr.get_structure(path, name, force_num_chains); } )
            .def("add_motif",[] (Manager  & ptr, String const & path, String name, util::MotifType mtype) {
                ptr.add_motif(path, name, mtype); } )
            .def("add_motif",[] (Manager  & ptr, motif::MotifOP const & m, String name) {
                ptr.add_motif(m, name); } )
            .def("register_motif",[] (Manager  & ptr, motif::MotifOP const & m) {
                ptr.register_motif(m); } )
            .def("register_extra_motif_ensembles",[] (Manager  & ptr, String const & f_name) {
                ptr.register_extra_motif_ensembles(f_name); } )
            .def("has_supplied_motif_ensemble",[] (Manager  & ptr, String const & m_name, String const & end_name) -> int {
                return ptr.has_supplied_motif_ensemble(m_name, end_name); } )
            .def("get_supplied_motif_ensemble",[] (Manager  & ptr, String const & m_name, String const & end_name) -> motif::MotifEnsembleOP const &  {
                return ptr.get_supplied_motif_ensemble(m_name, end_name); } )
            .def("get_helix_name",[] (Manager  & ptr, String const & str) -> String {
                return ptr.get_helix_name(str); } )
            .def("get_motif_from_state",[] (Manager  & ptr, motif::MotifStateOP ms) -> motif::MotifOP {
                return ptr.get_motif_from_state(ms); } )
            ;

    py::class_<MotifEnsembleSqliteConnection, std::shared_ptr<MotifEnsembleSqliteConnection>>(m, "MotifEnsembleSqliteConnection")
            // ctors
            .def(py::init<>())
            .def(py::init<String const &>())
                    // methods
            .def("next",[] (MotifEnsembleSqliteConnection  & ptr) -> MotifEnsembleSqliteDataOP const & {
                return ptr.next(); } )
            ;

    py::class_<MotifEnsembleSqliteData, std::shared_ptr<MotifEnsembleSqliteData>>(m, "MotifEnsembleSqliteData")
            // ctors
            .def(py::init<>())
                    // public attributes
            .def_readwrite("data", &MotifEnsembleSqliteData::data)
            .def_readwrite("name", &MotifEnsembleSqliteData::name)
            .def_readwrite("id", &MotifEnsembleSqliteData::id)
            ;

    py::class_<MotifSqliteConnection, std::shared_ptr<MotifSqliteConnection>>(m, "MotifSqliteConnection")
            // ctors
            .def(py::init<>())
            .def(py::init<String const &>())
                    // methods
            .def("next",[] (MotifSqliteConnection  & ptr) -> MotifSqliteDataOP const & {
                return ptr.next(); } )
            .def("contains",[] (MotifSqliteConnection  & ptr) -> MotifSqliteDataOP const & {
                return ptr.contains(); } )
            .def("clear",[] (MotifSqliteConnection  & ptr) {
                ptr.clear(); } )
            ;

    py::class_<MotifSqliteData, std::shared_ptr<MotifSqliteData>>(m, "MotifSqliteData")
            // ctors
            .def(py::init<>())
                    // public attributes
            .def_readwrite("data", &MotifSqliteData::data)
            .def_readwrite("name", &MotifSqliteData::name)
            .def_readwrite("end_name", &MotifSqliteData::end_name)
            .def_readwrite("end_id", &MotifSqliteData::end_id)
            .def_readwrite("id", &MotifSqliteData::id)
            ;

    py::class_<MotifSqliteLibrary, std::shared_ptr<MotifSqliteLibrary>>(m, "MotifSqliteLibrary")
            // ctors
            .def(py::init<String const &>())
            .def(py::init<int,std::filesystem::path const &>())
                    // methods
            .def("begin",[] (MotifSqliteLibrary  & ptr) -> MotifSqliteLibrary::iterator {
                return ptr.begin(); } )
            .def("end",[] (MotifSqliteLibrary  & ptr) -> MotifSqliteLibrary::iterator {
                return ptr.end(); } )
            .def("get_libnames",[] (MotifSqliteLibrary  & ptr) -> StringStringMap {
                return ptr.get_libnames(); } )
            .def("get",[] (MotifSqliteLibrary  & ptr, String const & name, String const & end_id, String const & end_name, String const & id) -> motif::MotifOP {
                return ptr.get(name, end_id, end_name, id); } )
            .def("get_multi",[] (MotifSqliteLibrary  & ptr, String const & name, String const & end_id, String const & end_name, String const & id) -> motif::MotifOPs {
                return ptr.get_multi(name, end_id, end_name, id); } )
            .def("contains",[] (MotifSqliteLibrary  & ptr, String const & name, String const & end_id, String const & end_name, String const & id) -> int {
                return ptr.contains(name, end_id, end_name, id); } )
            .def("get_random",[] (MotifSqliteLibrary  & ptr) -> motif::MotifOP {
                return ptr.get_random(); } )
            .def("load_all",[] (MotifSqliteLibrary  & ptr, int limit) {
                ptr.load_all(limit); } )
            ;

    py::class_<MotifStateEnsembleSqliteLibrary, std::shared_ptr<MotifStateEnsembleSqliteLibrary>>(m, "MotifStateEnsembleSqliteLibrary")
            // ctors
            .def(py::init<String const &>())
                    // methods
            .def("begin",[] (MotifStateEnsembleSqliteLibrary  & ptr) -> MotifStateEnsembleSqliteLibrary::iterator {
                return ptr.begin(); } )
            .def("end",[] (MotifStateEnsembleSqliteLibrary  & ptr) -> MotifStateEnsembleSqliteLibrary::iterator {
                return ptr.end(); } )
            .def("get_libnames",[] (MotifStateEnsembleSqliteLibrary  & ptr) ->  StringStringMap {
                return ptr.get_libnames(); } )
            .def("get",[] (MotifStateEnsembleSqliteLibrary  & ptr, String const & name, String const & id) -> motif::MotifStateEnsembleOP {
                return ptr.get(name, id); } )
            .def("get_multi",[] (MotifStateEnsembleSqliteLibrary  & ptr, String const & name, String const & id) -> motif::MotifStateEnsembleOPs {
                return ptr.get_multi(name, id); } )
            .def("contains",[] (MotifStateEnsembleSqliteLibrary  & ptr, String const & name, String const & id) -> int {
                return ptr.contains(name, id); } )
            .def("get_random",[] (MotifStateEnsembleSqliteLibrary  & ptr) -> motif::MotifStateEnsembleOP {
                return ptr.get_random(); } )
            .def("load_all",[] (MotifStateEnsembleSqliteLibrary  & ptr, int limit) {
                ptr.load_all(limit); } )
            ;

    py::class_<MotifStateSqliteLibrary, std::shared_ptr<MotifStateSqliteLibrary>>(m, "MotifStateSqliteLibrary")
            // ctors
            .def(py::init<String const &>())
                    // methods
            .def("begin",[] (MotifStateSqliteLibrary  & ptr) -> MotifStateSqliteLibrary::iterator {
                return ptr.begin(); } )
            .def("end",[] (MotifStateSqliteLibrary  & ptr) -> MotifStateSqliteLibrary::iterator {
                return ptr.end(); } )
            .def("get_libnames",[] (MotifStateSqliteLibrary  & ptr) -> StringStringMap {
                return ptr.get_libnames(); } )
            .def("get",[] (MotifStateSqliteLibrary  & ptr, String const & name, String const & end_id, String const & end_name, String const & id) -> motif::MotifStateOP {
                return ptr.get(name, end_id, end_name, id); } )
            .def("get_multi",[] (MotifStateSqliteLibrary  & ptr, String const & name, String const & end_id, String const & end_name, String const & id) -> motif::MotifStateOPs {
                return ptr.get_multi(name, end_id, end_name, id); } )
            .def("contains",[] (MotifStateSqliteLibrary  & ptr, String const & name, String const & end_id, String const & end_name, String const & id) -> int {
                return ptr.contains(name, end_id, end_name, id); } )
            .def("get_random",[] (MotifStateSqliteLibrary  & ptr) -> motif::MotifStateOP {
                return ptr.get_random(); } )
            .def("load_all",[] (MotifStateSqliteLibrary  & ptr, int limit) {
                ptr.load_all(limit); } )
            .def("get_name",[] (MotifStateSqliteLibrary  & ptr) -> String const & {
                return ptr.get_name(); } )
            ;

    py::class_<SqliteLibrary, std::shared_ptr<SqliteLibrary>>(m, "SqliteLibrary")
            // ctors
            .def(py::init<>())
            ;
/*
        py::class_<iterator, std::shared_ptr<iterator>>(m, "iterator")
		// ctors
		.def(py::init<std::map<String, motif::MotifStateOP>::iterator const &>())
		// operators
		.def(py::self == py::self)
		.def(py::self != py::self)
		;

        py::class_<iterator, std::shared_ptr<iterator>>(m, "iterator")
		// ctors
		.def(py::init<std::map<String, motif::MotifStateEnsembleOP>::iterator const &>())
		// operators
		.def(py::self == py::self)
		.def(py::self != py::self)
		;

        py::class_<iterator, std::shared_ptr<iterator>>(m, "iterator")
		// ctors
		.def(py::init<std::map<String, motif::MotifOP>::iterator const &>())
		// operators
		.def(py::self == py::self)
		.def(py::self != py::self)
		;
*/

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// motif_data_structures
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // free functions

    m.def("graph_to_tree", [] (MotifGraphOP const & mg, data_structure::graph::GraphNodeOP<motif::MotifOP> start, structure::BasepairOP last_end) -> MotifTreeOP {
              return graph_to_tree(mg, start, last_end); },
          py::arg("mg"),
          py::arg("start") = nullptr,
          py::arg("last_end") = nullptr
    );

    m.def("motif_state_ensemble_graph_from_motif_graph", [] (MotifGraph & mg, resources::Manager & man, MotifStateEnsembleGraph & meg, std::map<int, int> & map) {
              motif_state_ensemble_graph_from_motif_graph(mg, man, meg, map); },
          py::arg("graph"),
          py::arg("manager"),
          py::arg("ensemble_graph"),
          py::arg("map")
    );
    // classes

    py::class_<ChainNodeData, std::shared_ptr<ChainNodeData>>(m, "ChainNodeData")
            // ctors
            .def(py::init<>())
            .def(py::init<structure::ChainOP const &,util::Uuid const &,int,int>())
            .def(py::init<ChainNodeData const &>())
                    // methods
            .def("included_res",[] (ChainNodeData  & ptr) -> structure::ResidueOPs {
                return ptr.included_res(); } )
                    // public attributes
            .def_readwrite("c", &ChainNodeData::c)
            .def_readwrite("m_id", &ChainNodeData::m_id)
            .def_readwrite("prime5_override", &ChainNodeData::prime5_override)
            .def_readwrite("prime3_override", &ChainNodeData::prime3_override)
            ;

    py::class_<GraphtoTree, std::shared_ptr<GraphtoTree>>(m, "GraphtoTree")
            // ctors
            .def(py::init<>())
                    // methods
            .def("convert",[] (GraphtoTree  & ptr, MotifGraphOP const & mg, data_structure::graph::GraphNodeOP<motif::MotifOP> start, int start_end_index, data_structure::graph::GraphNodeOP<motif::MotifOP> last_node) -> MotifTreeOP {
                return ptr.convert(mg, start, start_end_index, last_node); } )
            ;

    py::class_<MSNodeData, std::shared_ptr<MSNodeData>>(m, "MSNodeData")
            // ctors
            .def(py::init<motif::MotifStateOP const &>())
            .def(py::init<MSNodeData const &>())
                    // methods
            .def("get_end_state",[] (MSNodeData  & ptr, String const & name) ->  structure::BasepairStateOP {
                return ptr.get_end_state(name); } )
            .def("get_end_state",[] (MSNodeData  & ptr, int i) ->  structure::BasepairStateOP {
                return ptr.get_end_state(i); } )
            .def("get_end_index",[] (MSNodeData  & ptr, String const & name) ->  int {
                return ptr.get_end_index(name); } )
            .def("name",[] (MSNodeData  & ptr) ->  String const & {
                return ptr.name(); } )
            .def("block_end_add",[] (MSNodeData  & ptr) ->  int {
                return ptr.block_end_add(); } )
            .def("end_name",[] (MSNodeData  & ptr, int i) ->  String const & {
                return ptr.end_name(i); } )
            .def("uuid",[] (MSNodeData  & ptr) ->  util::Uuid const & {
                return ptr.uuid(); } )
            .def("uuid",[] (MSNodeData  & ptr, util::Uuid const & uuid) {
                ptr.uuid(uuid); } )
                    // public attributes
            .def_readwrite("ref_state", &MSNodeData::ref_state)
            .def_readwrite("cur_state", &MSNodeData::cur_state)
            ;

    py::class_<MotifConnection, std::shared_ptr<MotifConnection>>(m, "MotifConnection")
            // ctors
            .def(py::init<>())
            .def(py::init<int,int,String const &,String const &>())
            .def(py::init<MotifConnection const &>())
                    // methods
            .def("i",[] (MotifConnection  & ptr) ->  int {
                return ptr.i(); } )
            .def("j",[] (MotifConnection  & ptr) ->  int {
                return ptr.j(); } )
            .def("name_i",[] (MotifConnection  & ptr) ->  String const & {
                return ptr.name_i(); } )
            .def("name_j",[] (MotifConnection  & ptr) ->  String const & {
                return ptr.name_j(); } )
            .def("name_i",[] (MotifConnection  & ptr, String const & name) {
                ptr.name_i(name); } )
            .def("name_j",[] (MotifConnection  & ptr, String const & name) {
                ptr.name_j(name); } )
            .def("to_str",[] (MotifConnection  & ptr) -> String {
                return ptr.to_str(); } )
            ;

    py::class_<MotifConnections, std::shared_ptr<MotifConnections>>(m, "MotifConnections")
            // ctors
            .def(py::init<>())
            .def(py::init<MotifConnections const>())
                    // methods
            .def("begin",[] (MotifConnections  & ptr) -> MotifConnectionOPs::iterator {
                return ptr.begin(); } )
            .def("end",[] (MotifConnections  & ptr) -> MotifConnectionOPs::iterator {
                return ptr.end(); } )
            .def("begin",[] (MotifConnections const & ptr) -> MotifConnectionOPs::const_iterator {
                return ptr.begin(); } )
            .def("end",[] (MotifConnections const & ptr) -> MotifConnectionOPs::const_iterator {
                return ptr.end(); } )
            .def("size",[] (MotifConnections const & ptr) ->  size_t {
                return ptr.size(); } )
            .def("add_connection",[] (MotifConnections  & ptr, int i, int j, String const & name_i, String const & name_j) {
                ptr.add_connection(i, j, name_i, name_j); } )
            .def("remove_connections_to",[] (MotifConnections  & ptr, int index) {
                ptr.remove_connections_to(index); } )
            .def("in_connection",[] (MotifConnections const & ptr, int index, String const & name) -> bool {
                return ptr.in_connection(index, name); } )
            .def("update_connection_name",[] (MotifConnections  & ptr, int index, String const & name, String const & new_name) {
                ptr.update_connection_name(index, name, new_name); } )
            ;

    py::class_<MotifGraph, std::shared_ptr<MotifGraph>>(m, "MotifGraph")
            // ctors
            .def(py::init<>())
                    //.def(py::init<String const>())
            .def(py::init<String const,MotifGraphStringType const>())
            .def(py::init<MotifGraph const>())
                    // methods
            .def("update_indexes",[] (MotifGraph  & ptr, std::map<int, int> const & index_hash) {
                ptr.update_indexes(index_hash); } )
            .def("begin",[] (MotifGraph  & ptr) -> MotifGraph::iterator {
                return ptr.begin(); } )
            .def("end",[] (MotifGraph  & ptr) -> MotifGraph::iterator {
                return ptr.end(); } )
            .def("add_motif",[] (MotifGraph  & ptr, motif::MotifOP const & m, int parent_index, int parent_end_index, int orphan) -> int {
                return ptr.add_motif(m, parent_index, parent_end_index, orphan); } )
            .def("add_motif",[] (MotifGraph  & ptr, motif::MotifOP const & m, int parent_index, String const & p_end_name) -> int {
                return ptr.add_motif(m, parent_index, p_end_name); } )
            .def("add_motif_tree",[] (MotifGraph  & ptr, MotifTreeOP const & mt, int parent_index, int parent_end_index) {
                ptr.add_motif_tree(mt, parent_index, parent_end_index); } )
            .def("add_motif_tree",[] (MotifGraph  & ptr, MotifTreeOP const & mt, int parent_index, String const & parent_end_name) {
                ptr.add_motif_tree(mt, parent_index, parent_end_name); } )
            .def("add_motif_graph",[] (MotifGraph  & ptr, MotifGraph & mg, int parent_index, int parent_end_index) {
                ptr.add_motif_graph(mg, parent_index, parent_end_index); } )
            .def("add_connection",[] (MotifGraph  & ptr, int i, int j, String const & i_bp_name, String const & j_bp_name) {
                ptr.add_connection(i, j, i_bp_name, j_bp_name); } )
            .def("parent_index",[] (MotifGraph  & ptr, int node_index) ->  int {
                return ptr.parent_index(node_index); } )
            .def("parent_end_index",[] (MotifGraph  & ptr, int node_index) ->  int {
                return ptr.parent_end_index(node_index); } )
            .def("remove_motif",[] (MotifGraph  & ptr, int index) {
                ptr.remove_motif(index); } )
            .def("remove_level",[] (MotifGraph  & ptr, int level) {
                ptr.remove_level(level); } )
            .def("size",[] (MotifGraph  & ptr) ->  size_t {
                return ptr.size(); } )
            .def("increase_level",[] (MotifGraph  & ptr) {
                ptr.increase_level(); } )
            .def("last_node",[] (MotifGraph  & ptr) ->  data_structure::graph::GraphNodeOP<motif::MotifOP> const & {
                return ptr.last_node(); } )
            .def("oldest_node",[] (MotifGraph  & ptr) ->  data_structure::graph::GraphNodeOP<motif::MotifOP> {
                return ptr.oldest_node(); } )
            .def("get_node",[] (MotifGraph const & ptr, int i) ->  data_structure::graph::GraphNodeOP<motif::MotifOP> const & {
                return ptr.get_node(i); } )
            .def("get_node",[] (MotifGraph const & ptr, util::Uuid const & uuid) ->  data_structure::graph::GraphNodeOP<motif::MotifOP> const {
                return ptr.get_node(uuid); } )
            .def("get_node",[] (MotifGraph const & ptr, String const & m_name) ->  data_structure::graph::GraphNodeOP<motif::MotifOP> const {
                return ptr.get_node(m_name); } )
            .def("set_index",[] (MotifGraph  & ptr, int index) {
                ptr.set_index(index); } )
            .def("replace_ideal_helices",[] (MotifGraph  & ptr) {
                ptr.replace_ideal_helices(); } )
            .def("replace_ideal_helices",[] (MotifGraph  & ptr) {
                ptr.replace_ideal_helices(); } )
            .def("replace_helical_sequence",[] (MotifGraph  & ptr, secondary_structure::PoseOP const & ss) {
                ptr.replace_helical_sequence(ss); } )
            .def("replace_helical_sequence",[] (MotifGraph  & ptr, String const & seq) {
                ptr.replace_helical_sequence(seq); } )
            .def("designable_secondary_structure",[] (MotifGraph  & ptr) -> secondary_structure::PoseOP {
                return ptr.designable_secondary_structure(); } )
            .def("designable_sequence",[] (MotifGraph  & ptr) ->  String {
                return ptr.designable_sequence(); } )
            .def("write_pdbs",[] (MotifGraph  & ptr, String const & fname) {
                ptr.write_pdbs(fname); } )
            .def("topology_to_str",[] (MotifGraph  & ptr) -> String {
                return ptr.topology_to_str(); } )
            .def("to_str",[] (MotifGraph  & ptr) -> String {
                return ptr.to_str(); } )
            .def("_update_align_list",[] (MotifGraph  & ptr) {
                ptr._update_align_list(); } )
            .def("_update_merger",[] (MotifGraph  & ptr) {
                ptr._update_merger(); } )
            .def("beads",[] (MotifGraph  & ptr) ->  structure::Beads {
                return ptr.beads(); } )
            .def("get_available_end",[] (MotifGraph  & ptr, int end) -> structure::BasepairOP const & {
                return ptr.get_available_end(end); } )
            .def("get_available_end",[] (MotifGraph  & ptr, int pos) -> structure::BasepairOP const & {
                return ptr.get_available_end(pos); } )
            .def("get_available_end",[] (MotifGraph  & ptr, int pos) -> structure::BasepairOP const &  {
                return ptr.get_available_end(pos); } )
            .def("unaligned_nodes",[] (MotifGraph const & ptr) -> data_structure::graph::GraphNodeOPs<motif::MotifOP> const {
                return ptr.unaligned_nodes(); } )
            .def("aligned",[] (MotifGraph  & ptr) -> std::map<int, int> {
                return ptr.aligned(); } )
            .def("connections",[] (MotifGraph  & ptr) -> data_structure::graph::GraphConnectionOPs<motif::MotifOP> const & {
                return ptr.connections(); } )
            .def("get_structure",[] (MotifGraph  & ptr) ->  structure::RNAStructureOP const & {
                return ptr.get_structure(); } )
            .def("secondary_structure",[] (MotifGraph  & ptr) -> secondary_structure::PoseOP {
                return ptr.secondary_structure(); } )
            .def("to_pdb",[] (MotifGraph  & ptr, String const fname, int renumber, int close_chains, int conect_statement) {
                ptr.to_pdb(fname, renumber, close_chains, conect_statement); } )
            .def("pdb_str",[] (MotifGraph  & ptr, int renumber, int close_chains, int conect_statement) -> String {
                return ptr.pdb_str(renumber, close_chains, conect_statement); } )
            .def("sequence",[] (MotifGraph  & ptr) ->  String {
                return ptr.sequence(); } )
            .def("dot_bracket",[] (MotifGraph  & ptr) ->  String {
                return ptr.dot_bracket(); } )
            .def("get_int_option",[] (MotifGraph  & ptr, String const & name) ->  float {
                return ptr.get_int_option(name); } )
            .def("get_float_option",[] (MotifGraph  & ptr, String const & name) ->  float {
                return ptr.get_float_option(name); } )
            .def("get_string_option",[] (MotifGraph  & ptr, String const & name) ->  String {
                return ptr.get_string_option(name); } )
            .def("get_bool_option",[] (MotifGraph  & ptr, String const & name) ->  bool {
                return ptr.get_bool_option(name); } )
            .def("set_option_value",[] (MotifGraph  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (MotifGraph  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (MotifGraph  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (MotifGraph  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )
            ;

    py::class_<MotifMerger, std::shared_ptr<MotifMerger>>(m, "MotifMerger")
            // ctors
            .def(py::init<>())
            .def(py::init<MotifMerger const &,motif::MotifOPs const &>())
                    // methods
            .def("add_motif",[] (MotifMerger  & ptr, motif::MotifOP const & m) {
                ptr.add_motif(m); } )
            .def("add_motif",[] (MotifMerger  & ptr, motif::MotifOP const & m) {
                ptr.add_motif(m); } )
            .def("get_structure",[] (MotifMerger  & ptr) -> structure::RNAStructureOP const & {
                return ptr.get_structure(); } )
            .def("remove_motif",[] (MotifMerger  & ptr, motif::MotifOP const & m) {
                ptr.remove_motif(m); } )
            .def("update_motif",[] (MotifMerger  & ptr, motif::MotifOP const & m) {
                ptr.update_motif(m); } )
            .def("connect_motifs",[] (MotifMerger  & ptr, motif::MotifOP const & m1, motif::MotifOP const & m2, structure::BasepairOP const & m1_end, structure::BasepairOP const & m2_end) {
                ptr.connect_motifs(m1, m2, m1_end, m2_end); } )
            .def("secondary_structure",[] (MotifMerger  & ptr) -> secondary_structure::PoseOP {
                return ptr.secondary_structure(); } )
            .def("get_residue",[] (MotifMerger  & ptr, util::Uuid const & uuid) ->  structure::ResidueOP {
                return ptr.get_residue(uuid); } )
            .def("get_basepair",[] (MotifMerger  & ptr, util::Uuid const & uuid) ->  structure::BasepairOP {
                return ptr.get_basepair(uuid); } )
            .def("to_pdb",[] (MotifMerger  & ptr, String const & fname, int renumber, int close_chains, int conect_statements) {
                ptr.to_pdb(fname, renumber, close_chains, conect_statements); } )
            ;

    py::class_<MotifStateEnsembleGraph, std::shared_ptr<MotifStateEnsembleGraph>>(m, "MotifStateEnsembleGraph")
            // ctors
            .def(py::init<>())
            .def(py::init<MotifStateEnsembleGraph const &>())
                    // methods
            .def("begin",[] (MotifStateEnsembleGraph  & ptr) -> MotifStateEnsembleGraph::iterator {
                return ptr.begin(); } )
            .def("end",[] (MotifStateEnsembleGraph  & ptr) -> MotifStateEnsembleGraph::iterator {
                return ptr.end(); } )
            .def("begin",[] (MotifStateEnsembleGraph const & ptr) -> MotifStateEnsembleGraph::const_iterator {
                return ptr.begin(); } )
            .def("end",[] (MotifStateEnsembleGraph const & ptr) -> MotifStateEnsembleGraph::const_iterator {
                return ptr.end(); } )
            .def("size",[] (MotifStateEnsembleGraph  & ptr) -> size_t {
                return ptr.size(); } )
            .def("add_ensemble",[] (MotifStateEnsembleGraph  & ptr, motif::MotifStateEnsemble const & ensemble) -> int {
                return ptr.add_ensemble(ensemble); } )
            .def("add_ensemble",[] (MotifStateEnsembleGraph  & ptr, motif::MotifStateEnsemble const & ensemble, data_structure::NodeIndexandEdge const & parent_nie) -> int {
                return ptr.add_ensemble(ensemble, parent_nie); } )
//		.def("add_connection",[] (MotifStateEnsembleGraph  & ptr) -> int {
//		 return ptr.add_connection(); } )
            .def("get_ensemble",[] (MotifStateEnsembleGraph  & ptr, Index ni) ->  motif::MotifStateEnsemble const & {
                return ptr.get_ensemble(ni); } )
            .def("has_parent",[] (MotifStateEnsembleGraph const & ptr, Index ni) ->  bool {
                return ptr.has_parent(ni); } )
            .def("get_parent_index",[] (MotifStateEnsembleGraph const & ptr, Index ni) ->  Index {
                return ptr.get_parent_index(ni); } )
            .def("get_parent_end_index",[] (MotifStateEnsembleGraph const & ptr, Index ni) ->  Index {
                return ptr.get_parent_end_index(ni); } )
            .def("get_ensemble_connections",[] (MotifStateEnsembleGraph const & ptr, Index ni) ->  std::vector<data_structure::Edge const *> const & {
                return ptr.get_ensemble_connections(ni); } )
            .def("are_ensembles_connected",[] (MotifStateEnsembleGraph  & ptr, Index n1, Index n2) ->  bool {
                return ptr.are_ensembles_connected(n1, n2); } )
            ;

    py::class_<MotifStateEnsembleOPGraph, std::shared_ptr<MotifStateEnsembleOPGraph>>(m, "MotifStateEnsembleOPGraph")
            // ctors
            .def(py::init<>())
                    // methods
            .def("begin",[] (MotifStateEnsembleOPGraph  & ptr) -> MotifStateEnsembleOPGraph::iterator {
                return ptr.begin(); } )
            .def("end",[] (MotifStateEnsembleOPGraph  & ptr) -> MotifStateEnsembleOPGraph::iterator {
                return ptr.end(); } )
            .def("begin",[] (MotifStateEnsembleOPGraph const & ptr) -> MotifStateEnsembleOPGraph::const_iterator {
                return ptr.begin(); } )
            .def("end",[] (MotifStateEnsembleOPGraph const & ptr) -> MotifStateEnsembleOPGraph::const_iterator {
                return ptr.end(); } )
            .def("size",[] (MotifStateEnsembleOPGraph  & ptr) -> size_t {
                return ptr.size(); } )
            .def("add_ensemble",[] (MotifStateEnsembleOPGraph  & ptr, motif::MotifStateEnsembleOP ensemble) -> int {
                return ptr.add_ensemble(ensemble); } )
            .def("add_ensemble",[] (MotifStateEnsembleOPGraph  & ptr, motif::MotifStateEnsembleOP ensemble, data_structure::NodeIndexandEdge const & parent_nie) -> int {
                return ptr.add_ensemble(ensemble, parent_nie); } )
            .def("get_ensemble",[] (MotifStateEnsembleOPGraph  & ptr, Index ni) ->  motif::MotifStateEnsembleOP {
                return ptr.get_ensemble(ni); } )
            .def("has_parent",[] (MotifStateEnsembleOPGraph const & ptr, Index ni) ->  bool {
                return ptr.has_parent(ni); } )
            .def("get_parent_index",[] (MotifStateEnsembleOPGraph const & ptr, Index ni) ->  Index {
                return ptr.get_parent_index(ni); } )
            .def("get_parent_end_index",[] (MotifStateEnsembleOPGraph const & ptr, Index ni) ->  Index {
                return ptr.get_parent_end_index(ni); } )
            .def("get_leafs",[] (MotifStateEnsembleOPGraph  & ptr) ->  std::vector<data_structure::NodeIndexandEdge> {
                return ptr.get_leafs(); } )
            ;

    py::class_<MotifStateEnsembleTree, std::shared_ptr<MotifStateEnsembleTree>>(m, "MotifStateEnsembleTree")
            // ctors
            .def(py::init<>())
            .def(py::init<MotifTreeOP const>())
            .def(py::init<MotifStateTreeOP const>())
                    // methods
            .def("begin",[] (MotifStateEnsembleTree  & ptr) -> MotifStateEnsembleTree::iterator {
                return ptr.begin(); } )
            .def("end",[] (MotifStateEnsembleTree  & ptr) -> MotifStateEnsembleTree::iterator {
                return ptr.end(); } )
            .def("begin",[] (MotifStateEnsembleTree const & ptr) -> MotifStateEnsembleTree::const_iterator {
                return ptr.begin(); } )
            .def("end",[] (MotifStateEnsembleTree const & ptr) -> MotifStateEnsembleTree::const_iterator {
                return ptr.end(); } )
            .def("add_ensemble",[] (MotifStateEnsembleTree  & ptr, motif::MotifStateEnsembleOP const & ensemble, int parent_index, int parent_end_index) -> int {
                return ptr.add_ensemble(ensemble, parent_index, parent_end_index); } )
            .def("to_mst",[] (MotifStateEnsembleTree  & ptr) -> MotifStateTreeOP {
                return ptr.to_mst(); } )
            .def("size",[] (MotifStateEnsembleTree  & ptr) -> size_t {
                return ptr.size(); } )
            .def("get_node",[] (MotifStateEnsembleTree  & ptr, int i) -> MotifStateEnsembleTreeNodeOP const & {
                return ptr.get_node(i); } )
            .def("last_node",[] (MotifStateEnsembleTree  & ptr) ->  MotifStateEnsembleTreeNodeOP const & {
                return ptr.last_node(); } )
            ;

    py::class_<MotifStateEnsembleTreeEnumerator, std::shared_ptr<MotifStateEnsembleTreeEnumerator>>(m, "MotifStateEnsembleTreeEnumerator")
            // ctors
            .def(py::init<motif_data_structure::MotifStateEnsembleTreeOP>())
                    // methods
            .def("record",[] (MotifStateEnsembleTreeEnumerator  & ptr, String fname) {
                ptr.record(fname); } )
                    // public attributes
            .def_readwrite("mtst_", &MotifStateEnsembleTreeEnumerator::mtst_)
            ;

    py::class_<MotifStateGraph, std::shared_ptr<MotifStateGraph>>(m, "MotifStateGraph")
            // ctors
            .def(py::init<>())
            .def(py::init<MotifGraphOP const>())
            .def(py::init<MotifStateGraph const>())
                    // methods
            .def("begin",[] (MotifStateGraph  & ptr) -> MotifStateGraph::iterator {
                return ptr.begin(); } )
            .def("end",[] (MotifStateGraph  & ptr) -> MotifStateGraph::iterator {
                return ptr.end(); } )
            .def("size",[] (MotifStateGraph  & ptr) -> size_t {
                return ptr.size(); } )
            .def("add_state",[] (MotifStateGraph  & ptr, motif::MotifStateOP const & state, int parent_index, int parent_end_index, int orphan) -> int {
                return ptr.add_state(state, parent_index, parent_end_index, orphan); } )
            .def("add_state",[] (MotifStateGraph  & ptr, motif::MotifStateOP const & state, int parent_index, String const & parent_end_name) -> int {
                return ptr.add_state(state, parent_index, parent_end_name); } )
            .def("add_connection",[] (MotifStateGraph  & ptr, int i, int j, String const & i_bp_name, String const & j_bp_name) {
                ptr.add_connection(i, j, i_bp_name, j_bp_name); } )
            .def("replace_state",[] (MotifStateGraph  & ptr, int i, motif::MotifStateOP const & new_state) {
                ptr.replace_state(i, new_state); } )
            .def("remove_state",[] (MotifStateGraph  & ptr, int pos) {
                ptr.remove_state(pos); } )
            .def("remove_level",[] (MotifStateGraph  & ptr, int level) {
                ptr.remove_level(level); } )
            .def("last_node",[] (MotifStateGraph  & ptr) ->  data_structure::graph::GraphNodeOP<MSNodeDataOP> {
                return ptr.last_node(); } )
            .def("get_node",[] (MotifStateGraph const & ptr, int i) ->  data_structure::graph::GraphNodeOP<MSNodeDataOP> const & {
                return ptr.get_node(i); } )
            .def("get_node",[] (MotifStateGraph const & ptr, util::Uuid const & uuid) ->  data_structure::graph::GraphNodeOP<MSNodeDataOP> const {
                return ptr.get_node(uuid); } )
            .def("get_node",[] (MotifStateGraph const & ptr, String const & m_name) ->  data_structure::graph::GraphNodeOP<MSNodeDataOP> const {
                return ptr.get_node(m_name); } )
            .def("increase_level",[] (MotifStateGraph  & ptr) {
                ptr.increase_level(); } )
            .def("to_motif_graph",[] (MotifStateGraph  & ptr) -> MotifGraphOP {
                return ptr.to_motif_graph(); } )
            .def("unaligned_nodes",[] (MotifStateGraph const & ptr) -> data_structure::graph::GraphNodeOPs<MSNodeDataOP> const {
                return ptr.unaligned_nodes(); } )
            .def("get_int_option",[] (MotifStateGraph  & ptr, String const & name) ->  float {
                return ptr.get_int_option(name); } )
            .def("get_float_option",[] (MotifStateGraph  & ptr, String const & name) ->  float {
                return ptr.get_float_option(name); } )
            .def("get_string_option",[] (MotifStateGraph  & ptr, String const & name) ->  String {
                return ptr.get_string_option(name); } )
            .def("get_bool_option",[] (MotifStateGraph  & ptr, String const & name) ->  bool {
                return ptr.get_bool_option(name); } )
            .def("set_option_value",[] (MotifGraph  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (MotifGraph  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (MotifGraph  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (MotifGraph  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )
            ;

    py::class_<MotifStateTree, std::shared_ptr<MotifStateTree>>(m, "MotifStateTree")
            // ctors
            .def(py::init<>())
            .def(py::init<MotifTreeOP const>())
            .def(py::init<MotifStateTree const>())
            .def(py::init<String const>())
                    // methods
            .def("begin",[] (MotifStateTree  & ptr) -> MotifStateTree::iterator {
                return ptr.begin(); } )
            .def("end",[] (MotifStateTree  & ptr) -> MotifStateTree::iterator {
                return ptr.end(); } )
            .def("begin",[] (MotifStateTree const & ptr) -> MotifStateTree::const_iterator {
                return ptr.begin(); } )
            .def("end",[] (MotifStateTree const & ptr) -> MotifStateTree::const_iterator {
                return ptr.end(); } )
            .def("add_state",[] (MotifStateTree  & ptr, motif::MotifStateOP const & state, int parent_index, int parent_end_index) -> int {
                return ptr.add_state(state, parent_index, parent_end_index); } )
            .def("add_state",[] (MotifStateTree  & ptr, motif::MotifStateOP const & state, int parent_index, String const & parent_end_name) -> int {
                return ptr.add_state(state, parent_index, parent_end_name); } )
            .def("add_mst",[] (MotifStateTree  & ptr, MotifStateTreeOP const & mst, int parent_index, int parent_end_index) -> int {
                return ptr.add_mst(mst, parent_index, parent_end_index); } )
            .def("add_mst",[] (MotifStateTree  & ptr, MotifStateTreeOP const & mst, int parent_index, int parent_end_index) -> int {
                return ptr.add_mst(mst, parent_index, parent_end_index); } )
            .def("add_connection",[] (MotifStateTree  & ptr, int i, int j, String const & i_bp_name, String const & j_bp_name) {
                ptr.add_connection(i, j, i_bp_name, j_bp_name); } )
            .def("replace_state",[] (MotifStateTree  & ptr, int i, motif::MotifStateOP const & new_state) {
                ptr.replace_state(i, new_state); } )
            .def("remove_node",[] (MotifStateTree  & ptr, int i) {
                ptr.remove_node(i); } )
            .def("remove_node_level",[] (MotifStateTree  & ptr, int level) {
                ptr.remove_node_level(level); } )
            .def("to_motif_tree",[] (MotifStateTree  & ptr) -> MotifTreeOP {
                return ptr.to_motif_tree(); } )
            .def("topology_to_str",[] (MotifStateTree  & ptr) -> String {
                return ptr.topology_to_str(); } )
            .def("centers",[] (MotifStateTree  & ptr) ->  math::Points {
                return ptr.centers(); } )
            .def("connections",[] (MotifStateTree  & ptr) ->  MotifConnections const & {
                return ptr.connections(); } )
            .def("write_pdbs",[] (MotifStateTree  & ptr, String const & fname) {
                ptr.write_pdbs(fname); } )
            .def("get_structure",[] (MotifStateTree  & ptr) ->  structure::RNAStructureOP {
                return ptr.get_structure(); } )
            .def("size",[] (MotifStateTree  & ptr) -> size_t {
                return ptr.size(); } )
            .def("last_node",[] (MotifStateTree  & ptr) -> MotifStateTreeNodeOP const & {
                return ptr.last_node(); } )
            .def("get_node",[] (MotifStateTree  & ptr, int i) -> MotifStateTreeNodeOP const & {
                return ptr.get_node(i); } )
            .def("get_node",[] (MotifStateTree  & ptr, util::Uuid const & uuid) ->  MotifStateTreeNodeOP const & {
                return ptr.get_node(uuid); } )
            .def("get_node",[] (MotifStateTree  & ptr, String const & m_name) ->  MotifStateTreeNodeOP {
                return ptr.get_node(m_name); } )
            .def("increase_level",[] (MotifStateTree  & ptr) {
                ptr.increase_level(); } )
            .def("decrease_level",[] (MotifStateTree  & ptr) {
                ptr.decrease_level(); } )
            .def("get_int_option",[] (MotifStateTree  & ptr, String const & name) ->  float {
                return ptr.get_int_option(name); } )
            .def("get_float_option",[] (MotifStateTree  & ptr, String const & name) ->  float {
                return ptr.get_float_option(name); } )
            .def("get_string_option",[] (MotifStateTree  & ptr, String const & name) ->  String {
                return ptr.get_string_option(name); } )
            .def("get_bool_option",[] (MotifStateTree  & ptr, String const & name) ->  bool {
                return ptr.get_bool_option(name); } )
            .def("set_option_value",[] (MotifGraph  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (MotifGraph  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (MotifGraph  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (MotifGraph  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )

            ;

    py::class_<MotifTree, std::shared_ptr<MotifTree>>(m, "MotifTree")
            // ctors
            .def(py::init<>())
            .def(py::init<String const>())
            .def(py::init<String const,MotifTreeStringType>())
            .def(py::init<MotifTree const>())
                    // methods
            .def("begin",[] (MotifTree  & ptr) -> MotifTree::iterator {
                return ptr.begin(); } )
            .def("end",[] (MotifTree  & ptr) -> MotifTree::iterator {
                return ptr.end(); } )
            .def("begin",[] (MotifTree const & ptr) -> MotifTree::const_iterator {
                return ptr.begin(); } )
            .def("end",[] (MotifTree const & ptr) -> MotifTree::const_iterator {
                return ptr.end(); } )
            .def("add_motif",[] (MotifTree  & ptr, motif::MotifOP const & m, int parent_index, int parent_end_index) -> int {
                return ptr.add_motif(m, parent_index, parent_end_index); } )
            .def("add_motif",[] (MotifTree  & ptr, motif::MotifOP const & m, int parent_index, String parent_end_name) -> int {
                return ptr.add_motif(m, parent_index, parent_end_name); } )
            .def("add_connection",[] (MotifTree  & ptr, int i, int j, String const & i_bp_name, String const & j_bp_name) {
                ptr.add_connection(i, j, i_bp_name, j_bp_name); } )
            .def("remove_node",[] (MotifTree  & ptr, int i) {
                ptr.remove_node(i); } )
            .def("remove_node_level",[] (MotifTree  & ptr, int level) {
                ptr.remove_node_level(level); } )
            .def("to_pdb",[] (MotifTree  & ptr, String const fname, int renumber, int close_chains, int conect_statements) {
                ptr.to_pdb(fname, renumber, close_chains, conect_statements); } )
            .def("topology_to_str",[] (MotifTree  & ptr) -> String {
                return ptr.topology_to_str(); } )
            .def("to_pretty_str",[] (MotifTree  & ptr) -> String {
                return ptr.to_pretty_str(); } )
            .def("to_str",[] (MotifTree  & ptr) -> String {
                return ptr.to_str(); } )
            .def("_update_merger",[] (MotifTree  & ptr) {
                ptr._update_merger(); } )
            .def("connections",[] (MotifTree  & ptr) -> MotifConnections const & {
                return ptr.connections(); } )
            .def("beads",[] (MotifTree  & ptr) -> structure::Beads {
                return ptr.beads(); } )
            .def("size",[] (MotifTree  & ptr) -> size_t {
                return ptr.size(); } )
            .def("get_node",[] (MotifTree  & ptr, int i) ->  data_structure::tree::TreeNodeOP<motif::MotifOP> const & {
                return ptr.get_node(i); } )
            .def("get_node",[] (MotifTree  & ptr, util::Uuid const & uuid) ->  data_structure::tree::TreeNodeOP<motif::MotifOP> const & {
                return ptr.get_node(uuid); } )
            .def("get_node",[] (MotifTree  & ptr, String const & m_name) ->  data_structure::tree::TreeNodeOP<motif::MotifOP> {
                return ptr.get_node(m_name); } )
            .def("last_node",[] (MotifTree  & ptr) ->  data_structure::tree::TreeNodeOP<motif::MotifOP> const & {
                return ptr.last_node(); } )
            .def("write_pdbs",[] (MotifTree  & ptr, String const & fname) {
                ptr.write_pdbs(fname); } )
            .def("increase_level",[] (MotifTree  & ptr) {
                ptr.increase_level(); } )
            .def("decrease_level",[] (MotifTree  & ptr) {
                ptr.decrease_level(); } )
            .def("get_structure",[] (MotifTree  & ptr) ->  structure::RNAStructureOP const & {
                return ptr.get_structure(); } )
            .def("secondary_structure",[] (MotifTree  & ptr) -> secondary_structure::PoseOP {
                return ptr.secondary_structure(); } )
            .def("designable_secondary_structure",[] (MotifTree  & ptr) -> secondary_structure::PoseOP {
                return ptr.designable_secondary_structure(); } )
            .def("get_int_option",[] (MotifTree  & ptr, String const & name) ->  float {
                return ptr.get_int_option(name); } )
            .def("get_float_option",[] (MotifTree  & ptr, String const & name) ->  float {
                return ptr.get_float_option(name); } )
            .def("get_string_option",[] (MotifTree  & ptr, String const & name) ->  String {
                return ptr.get_string_option(name); } )
            .def("get_bool_option",[] (MotifTree  & ptr, String const & name) ->  bool {
                return ptr.get_bool_option(name); } )
            .def("set_option_value",[] (MotifGraph  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (MotifGraph  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (MotifGraph  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (MotifGraph  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )
            ;
/*
        py::class_<MotifTreePrinter, std::shared_ptr<MotifTreePrinter>>(m, "MotifTreePrinter")
		// methods
		.def("print_tree",[] (MotifTreePrinter  & ptr, MotifTree const & mt) -> String {
		 return ptr.print_tree(mt); } )
		;
        py::class_<_GraphtoTreeNode, std::shared_ptr<_GraphtoTreeNode>>(m, "_GraphtoTreeNode")
		// public attributes
		.def_readwrite("parent", &_GraphtoTreeNode::parent)
		.def_readwrite("parent_end_index", &_GraphtoTreeNode::parent_end_index)
		.def_readwrite("node", &_GraphtoTreeNode::node)
		.def_readwrite("motif", &_GraphtoTreeNode::motif)
		;

*/
    // enums
    py::enum_<MotifGraphStringType>(m,"MotifGraphStringType")
            .value("OLD",MotifGraphStringType::OLD)
            .value("MG",MotifGraphStringType::MG)
            .value("TOP",MotifGraphStringType::TOP)
            ;
    py::enum_<MotifTreeStringType>(m,"MotifTreeStringType")
            .value("MT_STR",MotifTreeStringType::MT_STR)
            .value("TOP_STR",MotifTreeStringType::TOP_STR)
            ;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// thermo_fluctuation
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// classes

        py::class_<FrameScorer, std::shared_ptr<FrameScorer>>(m, "FrameScorer")
		// ctors
		.def(py::init<>())
		// methods
		.def("score",[] (FrameScorer  & ptr, structure::BasepairStateOP & state_1, structure::BasepairStateOP & state_2) ->  float {
		 return ptr.score(state_1, state_2); } )
		// inherited methods
		.def("score",[] (ThermoFlucScorer  & ptr, structure::BasepairStateOP & state_1, structure::BasepairStateOP & state_2) ->  float {
		 return ptr.score(state_1, state_2); } )
		;

        py::class_<FrameScorerDevel, std::shared_ptr<FrameScorerDevel>>(m, "FrameScorerDevel")
		// ctors
		.def(py::init<>())
		// methods
		.def("score",[] (FrameScorerDevel  & ptr, structure::BasepairStateOP & state_1, structure::BasepairStateOP & state_2) ->  float {
		 return ptr.score(state_1, state_2); } )
		.def("weight_d",[] (FrameScorerDevel  & ptr, float weight_d) {
		ptr.weight_d(weight_d); } )
		.def("weight_r",[] (FrameScorerDevel  & ptr, float weight_r) {
		ptr.weight_r(weight_r); } )
		.def("weight_d",[] (FrameScorerDevel  & ptr) ->  float {
		 return ptr.weight_d(); } )
		.def("weight_r",[] (FrameScorerDevel  & ptr) ->  float {
		 return ptr.weight_r(); } )
		// inherited methods
		.def("score",[] (ThermoFlucScorer  & ptr, structure::BasepairStateOP & state_1, structure::BasepairStateOP & state_2) ->  float {
		 return ptr.score(state_1, state_2); } )
		;

        py::class_<RunningAverage, std::shared_ptr<RunningAverage>>(m, "RunningAverage")
		// ctors
		.def(py::init<>())
		// methods
		.def("Update",[] (RunningAverage  & ptr, double valIn) -> double {
		 return ptr.Update(valIn); } )
		.def("Get",[] (RunningAverage  & ptr) -> double {
		 return ptr.Get(); } )
		.def("Count",[] (RunningAverage  & ptr) -> size_t {
		 return ptr.Count(); } )
		.def("Reset",[] (RunningAverage  & ptr) {
		ptr.Reset(); } )
		;

        py::class_<SixDScorer, std::shared_ptr<SixDScorer>>(m, "SixDScorer")
		// ctors
		.def(py::init<String const &,structure::BasepairOP>())
		// methods
		.def("score",[] (SixDScorer  & ptr, structure::BasepairStateOP & state_1, structure::BasepairStateOP & state_2) ->  float {
		 return ptr.score(state_1, state_2); } )
		// inherited methods
		.def("score",[] (ThermoFlucScorer  & ptr, structure::BasepairStateOP & state_1, structure::BasepairStateOP & state_2) ->  float {
		 return ptr.score(state_1, state_2); } )
		;

        py::class_<ThermoFlucRelax, std::shared_ptr<ThermoFlucRelax>>(m, "ThermoFlucRelax")
		// ctors
		.def(py::init<>())
		// methods
//		.def("run_with_graph",[] (ThermoFlucRelax  & ptr) {
//		ptr.run_with_graph(); } )
		;

        py::class_<ThermoFlucSampler, std::shared_ptr<ThermoFlucSampler>>(m, "ThermoFlucSampler")
		// ctors
		.def(py::init<>())
		// methods
		.def("sample",[] (ThermoFlucSampler  & ptr, int steps) {
		ptr.sample(steps); } )
		.def("setup",[] (ThermoFlucSampler  & ptr, motif_data_structure::MotifStateEnsembleTreeOP const & mset) {
		ptr.setup(mset); } )
		.def("next",[] (ThermoFlucSampler  & ptr) -> int {
		 return ptr.next(); } )
//		.def("undo",[] (ThermoFlucSampler  & ptr) {
//		ptr.undo(); } )
		.def("to_pdb",[] (ThermoFlucSampler  & ptr, String fname, int renumber) {
		ptr.to_pdb(fname, renumber); } )
		.def("temperature",[] (ThermoFlucSampler  & ptr) ->  float {
		 return ptr.temperature(); } )
		.def("mst",[] (ThermoFlucSampler  & ptr) ->  motif_data_structure::MotifStateTreeOP {
		 return ptr.mst(); } )
		.def("temperature",[] (ThermoFlucSampler  & ptr, float const & temp) {
		ptr.temperature(temp); } )
		.def("randomized_start",[] (ThermoFlucSampler  & ptr, bool random_start) {
		ptr.randomized_start(random_start); } )
		;

        py::class_<ThermoFlucScorer, std::shared_ptr<ThermoFlucScorer>>(m, "ThermoFlucScorer")
		// ctors
		.def(py::init<>())
		// methods
		.def("score",[] (ThermoFlucScorer  & ptr, structure::BasepairStateOP & state_1, structure::BasepairStateOP & state_2) ->  float {
		 return ptr.score(state_1, state_2); } )
		;

        py::class_<ThermoFlucSimulation, std::shared_ptr<ThermoFlucSimulation>>(m, "ThermoFlucSimulation")
		// ctors
		.def(py::init<>())
		// methods
		.def("setup",[] (ThermoFlucSimulation  & ptr, motif_data_structure::MotifStateEnsembleTreeOP const & mset, int ni1, int ni2, int ei1, int ei2) {
		ptr.setup(mset, ni1, ni2, ei1, ei2); } )
		.def("_check_sterics",[] (ThermoFlucSimulation  & ptr) ->  int {
		 return ptr._check_sterics(); } )
		.def("run",[] (ThermoFlucSimulation  & ptr) -> int {
		 return ptr.run(); } )
		.def("get_avg",[] (ThermoFlucSimulation  & ptr) -> double {
		 return ptr.get_avg(); } )
		.def("options",[] (ThermoFlucSimulation  & ptr) ->  base::Options & {
		 return ptr.options(); } )
		.def("get_int_option",[] (ThermoFlucSimulation  & ptr, String const & name) ->  float {
		 return ptr.get_int_option(name); } )
		.def("get_float_option",[] (ThermoFlucSimulation  & ptr, String const & name) ->  float {
		 return ptr.get_float_option(name); } )
		.def("get_string_option",[] (ThermoFlucSimulation  & ptr, String const & name) ->  String {
		 return ptr.get_string_option(name); } )
		.def("get_bool_option",[] (ThermoFlucSimulation  & ptr, String const & name) ->  bool {
		 return ptr.get_bool_option(name); } )
		.def("set_option_value",[] (ThermoFlucSimulation  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
		.def("set_option_value",[] (ThermoFlucSimulation  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
		.def("set_option_value",[] (ThermoFlucSimulation  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
		.def("set_option_value",[] (ThermoFlucSimulation  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )
		.def("update_var_options",[] (ThermoFlucSimulation  & ptr) {
		ptr.update_var_options(); } )
		;


    // classes

    py::class_<thermo_fluctuation::graph::FrameScorer, std::shared_ptr<thermo_fluctuation::graph::FrameScorer>>(m, "GraphFrameScorer")
            // ctors
            .def(py::init<>())
                    // methods
            .def("clone",[] (thermo_fluctuation::graph::FrameScorer const & ptr) -> thermo_fluctuation::graph::Scorer * {
                return ptr.clone(); } )
            .def("score",[] (thermo_fluctuation::graph::FrameScorer  & ptr, structure::BasepairState const & state_1, structure::BasepairState const & state_2) ->  float {
                return ptr.score(state_1, state_2); } )
                    // inherited methods
            .def("clone",[] (thermo_fluctuation::graph::Scorer const & ptr) -> thermo_fluctuation::graph::Scorer * {
                return ptr.clone(); } )
            .def("setup",[] (thermo_fluctuation::graph::Scorer  & ptr, bool target_an_aligned_end) {
                ptr.setup(target_an_aligned_end); } )
            .def("score",[] (thermo_fluctuation::graph::Scorer  & ptr, structure::BasepairState const & bpstate1, structure::BasepairState const & bpstate2 ) ->  float {
                return ptr.score(bpstate1, bpstate2); } )
            ;

    py::class_<thermo_fluctuation::graph::OldFrameScorer, std::shared_ptr<thermo_fluctuation::graph::OldFrameScorer>>(m, "OldFrameScorer")
            // ctors
            .def(py::init<>())
                    // methods
            .def("clone",[] (thermo_fluctuation::graph::OldFrameScorer const & ptr) -> thermo_fluctuation::graph::Scorer * {
                return ptr.clone(); } )
            .def("score",[] (thermo_fluctuation::graph::OldFrameScorer  & ptr, structure::BasepairState const & state_1, structure::BasepairState const & state_2) ->  float {
                return ptr.score(state_1, state_2); } )
                    // inherited methods
            .def("clone",[] (thermo_fluctuation::graph::Scorer const & ptr) -> thermo_fluctuation::graph::Scorer * {
                return ptr.clone(); } )
            .def("setup",[] (thermo_fluctuation::graph::Scorer  & ptr, bool target_an_aligned_end) {
                ptr.setup(target_an_aligned_end); } )
            .def("score",[] (thermo_fluctuation::graph::Scorer  & ptr, structure::BasepairState const & bpstate1, structure::BasepairState const & bpstate2) ->  float {
                return ptr.score(bpstate1, bpstate2); } )
            ;
/*
        py::class_<thermo_fluctuation::graph::Parameters, std::shared_ptr<thermo_fluctuation::graph::Parameters>>(m, "Parameters")
		// public attributes
		.def_readwrite("temperature", &Parameters::temperature)
		.def_readwrite("steric_radius", &Parameters::steric_radius)
		.def_readwrite("cutoff", &Parameters::cutoff)
		;
*/
    py::class_<thermo_fluctuation::graph::Sampler, std::shared_ptr<thermo_fluctuation::graph::Sampler>>(m, "Sampler")
            // ctors
            .def(py::init<motif_data_structure::MotifStateEnsembleGraph const &>())
                    // methods
            .def("get_initial_state",[] (thermo_fluctuation::graph::Sampler  & ptr) -> motif_data_structure::MotifStateGraphOP {
                return ptr.get_initial_state(); } )
            .def("next",[] (thermo_fluctuation::graph::Sampler  & ptr, motif_data_structure:: MotifStateGraphOP graph) -> int {
                return ptr.next(graph); } )
//		.def("undo",[] (Sampler  & ptr) {
//'		ptr.undo(); } )
            .def("set_temperature",[] (thermo_fluctuation::graph::Sampler  & ptr, float temp) {
                ptr.set_temperature(temp); } )
            ;
/*
        py::class_<Scorer, std::shared_ptr<Scorer>>(m, "Scorer")
		// ctors
		.def(py::init<>())
		// methods
		.def("clone",[] (Scorer const & ptr) -> Scorer * {
		 return ptr.clone(); } )
		.def("setup",[] (Scorer  & ptr, bool target_an_aligned_end) {
		ptr.setup(target_an_aligned_end); } )
		.def("score",[] (Scorer  & ptr, structure::BasepairState const &, structure::BasepairState const &) ->  float {
		 return ptr.score(&, &); } )
		;
*/
    py::class_<thermo_fluctuation::graph::Simulation, std::shared_ptr<thermo_fluctuation::graph::Simulation>>(m, "Simulation")
            // ctors
            .def(py::init<thermo_fluctuation::graph::ScorerOP,thermo_fluctuation::graph::sterics::StericsOP>())
                    // methods
            .def("setup",[] (thermo_fluctuation::graph::Simulation  & ptr, motif_data_structure::MotifStateEnsembleGraph const & mseg, data_structure::NodeIndexandEdge const & start, data_structure::NodeIndexandEdge const & end) {
                ptr.setup(mseg, start, end); } )
            .def("next",[] (thermo_fluctuation::graph::Simulation  & ptr) -> bool {
                return ptr.next(); } )
            .def("write_pdbs",[] (thermo_fluctuation::graph::Simulation  & ptr, String const & name) {
                ptr.write_pdbs(name); } )
            .def("get_pdb_str",[] (thermo_fluctuation::graph::Simulation  & ptr) -> String {
                return ptr.get_pdb_str(); } )
            .def("get_motif_graph",[] (thermo_fluctuation::graph::Simulation  & ptr) -> motif_data_structure::MotifGraphOP {
                return ptr.get_motif_graph(); } )
            .def("get_score",[] (thermo_fluctuation::graph::Simulation  & ptr) -> float {
                return ptr.get_score(); } )
            .def("get_int_option",[] (thermo_fluctuation::graph::Simulation  & ptr, String const & name) ->  float {
                return ptr.get_int_option(name); } )
            .def("get_float_option",[] (thermo_fluctuation::graph::Simulation  & ptr, String const & name) ->  float {
                return ptr.get_float_option(name); } )
            .def("get_string_option",[] (thermo_fluctuation::graph::Simulation  & ptr, String const & name) ->  String {
                return ptr.get_string_option(name); } )
            .def("get_bool_option",[] (thermo_fluctuation::graph::Simulation  & ptr, String const & name) ->  bool {
                return ptr.get_bool_option(name); } )
            .def("set_option_value",[] (thermo_fluctuation::graph::Simulation  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (thermo_fluctuation::graph::Simulation  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (thermo_fluctuation::graph::Simulation  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (thermo_fluctuation::graph::Simulation  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )
            ;

    py::class_<thermo_fluctuation::graph::logging::TargetBPInfoLogger, std::shared_ptr<thermo_fluctuation::graph::logging::TargetBPInfoLogger>>(m, "TargetBPInfoLogger")
            // ctors
            .def(py::init<String const &>())
                    // methods
            .def("clone",[] (thermo_fluctuation::graph::logging::TargetBPInfoLogger const & ptr) -> thermo_fluctuation::graph::logging::Logger * {
                return ptr.clone(); } )
            .def("setup",[] (thermo_fluctuation::graph::logging::TargetBPInfoLogger  & ptr, motif_data_structure::MotifStateGraphOP msg, data_structure::NodeIndexandEdge const & start, data_structure::NodeIndexandEdge const & end) {
                ptr.setup(msg, start, end); } )
            .def("log",[] (thermo_fluctuation::graph::logging::TargetBPInfoLogger  & ptr, motif_data_structure::MotifStateGraphOP msg, float score) {
                ptr.log(msg, score); } )
                    // inherited methods
            .def("clone",[] (thermo_fluctuation::graph::logging::Logger const & ptr) -> thermo_fluctuation::graph::logging::Logger * {
                return ptr.clone(); } )
            .def("setup",[] (thermo_fluctuation::graph::logging::Logger  & ptr, motif_data_structure:: MotifStateGraphOP graph, data_structure::NodeIndexandEdge const & nandi1, data_structure::NodeIndexandEdge const &nandi2) {
                ptr.setup(graph, nandi1, nandi2); } )
            .def("log",[] (thermo_fluctuation::graph::logging::Logger  & ptr, motif_data_structure:: MotifStateGraphOP graph, float  value) {
                ptr.log(graph, value); } )
            ;


    // classes

    py::class_<thermo_fluctuation::graph::sterics::NoSterics, std::shared_ptr<thermo_fluctuation::graph::sterics::NoSterics>>(m, "NoSterics")
            // ctors
            .def(py::init<>())
                    // methods
            .def("clone",[] (thermo_fluctuation::graph::sterics::NoSterics const & ptr) -> thermo_fluctuation::graph::sterics::Sterics * {
                return ptr.clone(); } )
            .def("clash",[] (thermo_fluctuation::graph::sterics::NoSterics  & ptr, motif_data_structure::MotifStateGraphOP msg) -> bool {
                return ptr.clash(msg); } )
                    // inherited methods
            .def("clone",[] (thermo_fluctuation::graph::sterics::Sterics const & ptr) -> thermo_fluctuation::graph::sterics::Sterics * {
                return ptr.clone(); } )
            .def("clash",[] (thermo_fluctuation::graph::sterics::Sterics  & ptr, motif_data_structure:: MotifStateGraphOP graph) -> bool {
                return ptr.clash(graph); } )
            ;

    py::class_<thermo_fluctuation::graph::sterics::SelectiveSterics, std::shared_ptr<thermo_fluctuation::graph::sterics::SelectiveSterics>>(m, "SelectiveSterics")
            // ctors
            .def(py::init<Ints const,Ints const,float>())
                    // methods
            .def("clone",[] (thermo_fluctuation::graph::sterics::SelectiveSterics const & ptr) -> thermo_fluctuation::graph::sterics::Sterics * {
                return ptr.clone(); } )
            .def("clash",[] (thermo_fluctuation::graph::sterics::SelectiveSterics  & ptr, motif_data_structure::MotifStateGraphOP msg) -> bool {
                return ptr.clash(msg); } )
                    // inherited methods
            .def("clone",[] (thermo_fluctuation::graph::sterics::Sterics const & ptr) -> thermo_fluctuation::graph::sterics::Sterics * {
                return ptr.clone(); } )
            .def("clash",[] (thermo_fluctuation::graph::sterics::Sterics  & ptr, motif_data_structure:: MotifStateGraphOP graph) -> bool {
                return ptr.clash(graph); } )
            ;
/*
        py::class_<Sterics, std::shared_ptr<Sterics>>(m, "Sterics")
		// ctors
		.def(py::init<>())
		// methods
		.def("clone",[] (Sterics const & ptr) -> Sterics * {
		 return ptr.clone(); } )
		.def("clash",[] (Sterics  & ptr, motif_data_structure:: MotifStateGraphOP graph) -> bool {
		 return ptr.clash(graph); } )
		;
*/
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// motif_search
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // classes

    py::class_<MotifStateMonteCarlo, std::shared_ptr<MotifStateMonteCarlo>>(m, "MotifStateMonteCarlo")
            // ctors
            .def(py::init<std::vector<motif::MotifStateOPs> const &>())
                    // methods
            .def("setup",[] (MotifStateMonteCarlo  & ptr, motif_data_structure::MotifStateGraphOP msg, int arg1, int arg2, int arg3, int arg4, bool arg5) {
                ptr.setup(msg, arg1, arg2, arg3, arg4, arg5); } )
            .def("run",[] (MotifStateMonteCarlo  & ptr) {
                ptr.run(); } )
            .def("start",[] (MotifStateMonteCarlo  & ptr) {
                ptr.start(); } )
            .def("next",[] (MotifStateMonteCarlo  & ptr) -> MotifStateMonteCarloSolutionOP {
                return ptr.next(); } )
            .def("next_state",[] (MotifStateMonteCarlo  & ptr) -> MotifStateMonteCarloSolutionNewOP {
                return ptr.next_state(); } )
            .def("finished",[] (MotifStateMonteCarlo  & ptr) -> bool {
                return ptr.finished(); } )
            .def("lookup",[] (MotifStateMonteCarlo  & ptr, util::StericLookupNew const & sl) {
                ptr.lookup(sl); } )
            .def("options",[] (MotifStateMonteCarlo  & ptr) ->  base::Options & {
                return ptr.options(); } )
            .def("get_int_option",[] (MotifStateMonteCarlo  & ptr, String const & name) ->  float {
                return ptr.get_int_option(name); } )
            .def("get_float_option",[] (MotifStateMonteCarlo  & ptr, String const & name) ->  float {
                return ptr.get_float_option(name); } )
            .def("get_string_option",[] (MotifStateMonteCarlo  & ptr, String const & name) ->  String {
                return ptr.get_string_option(name); } )
            .def("get_bool_option",[] (MotifStateMonteCarlo  & ptr, String const & name) ->  bool {
                return ptr.get_bool_option(name); } )
            .def("has_option",[] (MotifStateMonteCarlo  & ptr, String const & name) ->  bool {
                return ptr.has_option(name); } )
            .def("set_option_value",[] (MotifStateMonteCarlo  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (MotifStateMonteCarlo  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (MotifStateMonteCarlo  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (MotifStateMonteCarlo  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )
            ;

    py::class_<MotifStateMonteCarloSolution, std::shared_ptr<MotifStateMonteCarloSolution>>(m, "MotifStateMonteCarloSolution")
            // ctors
            .def(py::init<motif_data_structure::MotifGraphOP,float>())
                    // public attributes
            .def_readwrite("mg", &MotifStateMonteCarloSolution::mg)
            .def_readwrite("score", &MotifStateMonteCarloSolution::score)
            ;

    py::class_<MotifStateMonteCarloSolutionNew, std::shared_ptr<MotifStateMonteCarloSolutionNew>>(m, "MotifStateMonteCarloSolutionNew")
            // ctors
            .def(py::init<motif_data_structure::MotifStateGraphOP,float>())
                    // public attributes
            .def_readwrite("msg", &MotifStateMonteCarloSolutionNew::msg)
            .def_readwrite("score", &MotifStateMonteCarloSolutionNew::score)
            ;

    py::class_<NoExclusionFilter, std::shared_ptr<NoExclusionFilter>>(m, "NoExclusionFilter")
            // ctors
            .def(py::init<>())
                    // methods
            .def("clone",[] (NoExclusionFilter const & ptr) -> SolutionFilter * {
                return ptr.clone(); } )
            .def("accept",[] (NoExclusionFilter  & ptr, Strings const & motif_names) -> bool {
                return ptr.accept(motif_names); } )
                    // inherited methods
            .def("clone",[] (SolutionFilter const & ptr) -> SolutionFilter * {
                return ptr.clone(); } )
            .def("accept",[] (SolutionFilter  & ptr, Strings const & sols) -> bool {
                return ptr.accept(sols); } )
            ;

/*
        py::class_<Node, std::shared_ptr<Node>>(m, "Node")
		// methods
		.def("get_type",[] (Node const & ptr) ->  NodeType {
		 return ptr.get_type(); } )
		.def("get_lib_name",[] (Node const & ptr) ->  String const & {
		 return ptr.get_lib_name(); } )
		.def("get_motif_state",[] (Node const & ptr) ->  motif::MotifStateOP {
		 return ptr.get_motif_state(); } )
		.def("get_motif_state_ensemble",[] (Node const & ptr) ->  motif::MotifStateEnsembleOP {
		 return ptr.get_motif_state_ensemble(); } )
		;
        py::enum_<NodeType, std::shared_ptr<NodeType>>(m, "NodeType")

		;
        py::class_<path_finding::Parameters, std::shared_ptr<path_finding::Parameters>>(m, "Parameters")
		// public attributes
		.def_readwrite("max_helix_size", &Parameters::max_helix_size)
		.def_readwrite("min_helix_size", &path_finding::Parameters::min_helix_size)
		;
*/
    py::class_<motif_search::path_finding::Parameters, std::shared_ptr<motif_search::path_finding::Parameters>>(m, "motif_search::path_finding::Parameters")
            // public attributes
            .def_readwrite("sterics", &motif_search::path_finding::Parameters::sterics)
            .def_readwrite("helix_end", &motif_search::path_finding::Parameters::helix_end)
            .def_readwrite("max_node_level", &motif_search::path_finding::Parameters::max_node_level)
            .def_readwrite("min_size", &motif_search::path_finding::Parameters::min_size)
            .def_readwrite("max_size", &motif_search::path_finding::Parameters::max_size)
            .def_readwrite("max_solutions", &motif_search::path_finding::Parameters::max_solutions)
            .def_readwrite("min_node_level", &motif_search::path_finding::Parameters::min_node_level)
            .def_readwrite("accept_score", &motif_search::path_finding::Parameters::accept_score)
            .def_readwrite("min_ss_score", &motif_search::path_finding::Parameters::min_ss_score)
            .def_readwrite("max_steps", &motif_search::path_finding::Parameters::max_steps)
            .def_readwrite("return_best", &motif_search::path_finding::Parameters::return_best)
            ;
    py::class_<Problem, std::shared_ptr<Problem>>(m, "Problem")
            // ctors
            .def(py::init<structure::BasepairStateOP,structure::BasepairStateOP,util::StericLookupNewOP,bool>())
                    // public attributes
            .def_readwrite("start", &Problem::start)
            .def_readwrite("end", &Problem::end)
            .def_readwrite("lookup", &Problem::lookup)
            .def_readwrite("target_an_aligned_end", &Problem::target_an_aligned_end)
            ;

    py::class_<RemoveDuplicateHelices, std::shared_ptr<RemoveDuplicateHelices>>(m, "RemoveDuplicateHelices")
            // ctors
            .def(py::init<>())
                    // methods
            .def("clone",[] (RemoveDuplicateHelices const & ptr) -> SolutionFilter * {
                return ptr.clone(); } )
            .def("accept",[] (RemoveDuplicateHelices  & ptr, Strings const & motif_names) -> bool {
                return ptr.accept(motif_names); } )
                    // inherited methods
            .def("clone",[] (SolutionFilter const & ptr) -> SolutionFilter * {
                return ptr.clone(); } )
            .def("accept",[] (SolutionFilter  & ptr, Strings const & strs) -> bool {
                return ptr.accept(strs); } )
            ;
/*
        py::class_<motif_search::Search, std::shared_ptr<motif_search::Search>>(m, "motif_search::Search")
		// ctors
		.def(py::init<String const &>())
		// methods
		.def("clone",[] (motif_search::Search const & ptr) -> motif_search::Search * {
		 return ptr.clone(); } )
		.def("setup",[] (motif_search::Search  & ptr, motif_search::ProblemOP ) {
		ptr.setup(); } )
		.def("start",[] (motif_search::Search  & ptr) {
		ptr.start(); } )
		.def("finished",[] (motif_search::Search  & ptr) -> bool {
		 return ptr.finished(); } )
		.def("next",[] (motif_search::Search  & ptr) -> SolutionOP {
		 return ptr.next(); } )
		.def("name",[] (motif_search::Search  & ptr) -> String const & {
		 return ptr.name(); } )
		.def("get_int_option",[] (motif_search::Search  & ptr, String const & name) ->  float {
		 return ptr.get_int_option(name); } )
		.def("get_float_option",[] (motif_search::Search  & ptr, String const & name) ->  float {
		 return ptr.get_float_option(name); } )
		.def("get_string_option",[] (motif_search::Search  & ptr, String const & name) ->  String {
		 return ptr.get_string_option(name); } )
		.def("get_bool_option",[] (motif_search::Search  & ptr, String const & name) ->  bool {
		 return ptr.get_bool_option(name); } )
		.def("set_option_value",[] (motif_search::Search  & ptr, String const & name, T const & val) {
		ptr.set_option_value(name, val); } )
		;
*/
    py::class_<Solution, std::shared_ptr<Solution>>(m, "Solution")
            // ctors
            .def(py::init<motif_data_structure::MotifStateGraphOP,float>())
                    // public attributes
            .def_readwrite("graph", &Solution::graph)
            .def_readwrite("score", &Solution::score)
            ;
/*
        py::class_<SolutionFilter, std::shared_ptr<SolutionFilter>>(m, "SolutionFilter")
		// ctors
		.def(py::init<>())
		// methods
		.def("clone",[] (SolutionFilter const & ptr) -> SolutionFilter * {
		 return ptr.clone(); } )
		.def("accept",[] (SolutionFilter  & ptr, Strings const &) -> bool {
		 return ptr.accept(&); } )
		;
*/
    py::class_<SolutionToplogy, std::shared_ptr<SolutionToplogy>>(m, "SolutionToplogy")
            // ctors
            .def(py::init<motif_data_structure::MotifStateEnsembleOPGraphOP>())
                    // methods
            .def("begin",[] (SolutionToplogy  & ptr) -> SolutionToplogy::iterator {
                return ptr.begin(); } )
            .def("end",[] (SolutionToplogy  & ptr) -> SolutionToplogy::iterator {
                return ptr.end(); } )
            .def("begin",[] (SolutionToplogy const & ptr) -> SolutionToplogy::const_iterator {
                return ptr.begin(); } )
            .def("end",[] (SolutionToplogy const & ptr) -> SolutionToplogy::const_iterator {
                return ptr.end(); } )
            .def("initialize_solution",[] (SolutionToplogy  & ptr, structure::BasepairStateOP bp_state) -> motif_data_structure::MotifStateGraphOP {
                return ptr.initialize_solution(bp_state); } )
            .def("initialize_solution_no_start",[] (SolutionToplogy  & ptr, structure::BasepairStateOP bp_state) -> motif_data_structure::MotifStateGraphOP {
                return ptr.initialize_solution_no_start(bp_state); } )
            .def("get_motif_state",[] (SolutionToplogy  & ptr, Index pos) -> motif::MotifStateOP {
                return ptr.get_motif_state(pos); } )
            .def("get_solution_nie",[] (SolutionToplogy  & ptr) -> std::vector<data_structure::NodeIndexandEdge> const & {
                return ptr.get_solution_nie(); } )
            .def("size",[] (SolutionToplogy  & ptr) ->  size_t {
                return ptr.size(); } )
            .def("get_ensemble_size",[] (SolutionToplogy  & ptr, Index pos) ->  size_t {
                return ptr.get_ensemble_size(pos); } )
            ;

    py::class_<SolutionToplogyFactory, std::shared_ptr<SolutionToplogyFactory>>(m, "SolutionToplogyFactory")
            // ctors
            .def(py::init<>())
                    // methods
            .def("generate_toplogy",[] (SolutionToplogyFactory  & ptr, SolutionTopologyTemplate const & sol_template) -> SolutionToplogyOP {
                return ptr.generate_toplogy(sol_template); } )
            .def("get_int_option",[] (SolutionToplogyFactory  & ptr, String const & name) ->  float {
                return ptr.get_int_option(name); } )
            .def("get_float_option",[] (SolutionToplogyFactory  & ptr, String const & name) ->  float {
                return ptr.get_float_option(name); } )
            .def("get_string_option",[] (SolutionToplogyFactory  & ptr, String const & name) ->  String {
                return ptr.get_string_option(name); } )
            .def("get_bool_option",[] (SolutionToplogyFactory  & ptr, String const & name) ->  bool {
                return ptr.get_bool_option(name); } )
            .def("set_option_value",[] (MotifStateMonteCarlo  & ptr, String const & name, float const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (MotifStateMonteCarlo  & ptr, String const & name, int const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (MotifStateMonteCarlo  & ptr, String const & name, bool const & val) { ptr.set_option_value(name, val); } )
            .def("set_option_value",[] (MotifStateMonteCarlo  & ptr, String const & name, String const & val) { ptr.set_option_value(name, val); } )
            ;

    py::class_<SolutionTopologyTemplate, std::shared_ptr<SolutionTopologyTemplate>>(m, "SolutionTopologyTemplate")
            // ctors
            .def(py::init<>())
                    // methods
            .def("begin",[] (SolutionTopologyTemplate  & ptr) -> SolutionTopologyTemplate::iterator {
                return ptr.begin(); } )
            .def("end",[] (SolutionTopologyTemplate  & ptr) -> SolutionTopologyTemplate::iterator {
                return ptr.end(); } )
            .def("begin",[] (SolutionTopologyTemplate const & ptr) -> SolutionTopologyTemplate::const_iterator {
                return ptr.begin(); } )
            .def("end",[] (SolutionTopologyTemplate const & ptr) -> SolutionTopologyTemplate::const_iterator {
                return ptr.end(); } )
            .def("add_library",[] (SolutionTopologyTemplate  & ptr, String const & lib_name) {
                ptr.add_library(lib_name); } )
            .def("add_library",[] (SolutionTopologyTemplate  & ptr, String const & lib_name, data_structure::NodeIndexandEdge const & parent_nie) {
                ptr.add_library(lib_name, parent_nie); } )
            .def("add_motif_state",[] (SolutionTopologyTemplate  & ptr, motif::MotifStateOP ms) {
                ptr.add_motif_state(ms); } )
            .def("add_motif_state",[] (SolutionTopologyTemplate  & ptr, motif::MotifStateOP ms, data_structure::NodeIndexandEdge const & parent_nie) {
                ptr.add_motif_state(ms, parent_nie); } )
            .def("add_ensemble",[] (SolutionTopologyTemplate  & ptr, motif::MotifStateEnsembleOP mse) {
                ptr.add_ensemble(mse); } )
            .def("add_ensemble",[] (SolutionTopologyTemplate  & ptr, motif::MotifStateEnsembleOP mse, data_structure::NodeIndexandEdge const & parent_nie) {
                ptr.add_ensemble(mse, parent_nie); } )
            .def("has_parent",[] (SolutionTopologyTemplate const & ptr, Index ni) ->  bool {
                return ptr.has_parent(ni); } )
            .def("get_parent_index",[] (SolutionTopologyTemplate const & ptr, Index ni) ->  Index {
                return ptr.get_parent_index(ni); } )
            .def("get_parent_end_index",[] (SolutionTopologyTemplate const & ptr, Index ni) ->  Index {
                return ptr.get_parent_end_index(ni); } )
            ;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// motif_search::exhaustive
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // classes

    py::class_<motif_search::exhaustive::DefaultScorer, std::shared_ptr<motif_search::exhaustive::DefaultScorer>>(m, "ExhaustiveDefaultScorer")
            // ctors
            .def(py::init<>())
                    // methods
            .def("clone",[] (motif_search::exhaustive::DefaultScorer const & ptr) -> motif_search::exhaustive::Scorer * {
                return ptr.clone(); } )
            .def("score",[] (motif_search::exhaustive::DefaultScorer  & ptr, structure::BasepairState const & bps) ->  float {
                return ptr.score(bps); } )
                    // inherited methods
            .def("clone",[] (motif_search::exhaustive::Scorer const & ptr) -> motif_search::exhaustive::Scorer * {
                return ptr.clone(); } )
            .def("set_target",[] (motif_search::exhaustive::Scorer  & ptr, structure::BasepairStateOP target, bool target_an_aligned_end) {
                ptr.set_target(target, target_an_aligned_end); } )
            .def("score",[] (motif_search::exhaustive::Scorer  & ptr, structure::BasepairState const & bps) ->  float {
                return ptr.score(bps); } )
            ;

    py::class_<motif_search::exhaustive::MotifStateEnumerator, std::shared_ptr<motif_search::exhaustive::MotifStateEnumerator>>(m, "MotifStateEnumerator")
            // ctors
            .def(py::init<motif_search::SolutionToplogy>())
                    // methods
            .def("start",[] (motif_search::exhaustive::MotifStateEnumerator  & ptr, structure:: BasepairStateOP bp) {
                ptr.start(bp); } )
            .def("finished",[] (motif_search::exhaustive::MotifStateEnumerator  & ptr) -> bool {
                return ptr.finished(); } )
            .def("next",[] (motif_search::exhaustive::MotifStateEnumerator  & ptr) {
                ptr.next(); } )
            .def("top_state",[] (motif_search::exhaustive::MotifStateEnumerator  & ptr) -> motif::MotifStateOP {
                return ptr.top_state(); } )
            .def("all_states",[] (motif_search::exhaustive::MotifStateEnumerator  & ptr) -> motif::MotifStateOPs const & {
                return ptr.all_states(); } )
            .def("set_size_limit",[] (motif_search::exhaustive::MotifStateEnumerator  & ptr, int size_limit) {
                ptr.set_size_limit(size_limit); } )
            ;

/*
        py::class_<Scorer, std::shared_ptr<Scorer>>(m, "Scorer")
		// ctors
		.def(py::init<>())
		// methods
		.def("clone",[] (Scorer const & ptr) -> Scorer * {
		 return ptr.clone(); } )
		.def("set_target",[] (Scorer  & ptr, motif_search::exhaustivestructure::BasepairStateOP target, bool target_an_aligned_end) {
		ptr.set_target(target, target_an_aligned_end); } )
		.def("score",[] (Scorer  & ptr, structure::BasepairState const & bps) ->  float {
		 return ptr.score(bps); } )
		;
*/
    py::class_<motif_search::exhaustive::ScorerFactory, std::shared_ptr<motif_search::exhaustive::ScorerFactory>>(m, "ExhaustiveScorerFactory")
            // ctors
            .def(py::init<>())
                    // methods
            .def("get_scorer",[] (motif_search::exhaustive::ScorerFactory  & ptr, String const & scorer_name) -> motif_search::exhaustive::ScorerOP {
                return ptr.get_scorer(scorer_name); } )
            ;

    py::class_<motif_search::exhaustive::Search, std::shared_ptr<motif_search::exhaustive::Search>>(m, "ExhaustiveSearch")
            // ctors
            .def(py::init<motif_search::exhaustive::ScorerOP, motif_search::SolutionToplogy const & , motif_search::SolutionFilterOP>())
                    // methods
            .def("clone",[] (motif_search::exhaustive::Search const & ptr) -> motif_search::Search * {
                return ptr.clone(); } )
            .def("setup",[] (motif_search::exhaustive::Search  & ptr, motif_search::ProblemOP p) {
                ptr.setup(p); } )
            .def("start",[] (motif_search::exhaustive::Search  & ptr) {
                ptr.start(); } )
            .def("finished",[] (motif_search::exhaustive::Search  & ptr) -> bool {
                return ptr.finished(); } )
            .def("next",[] (motif_search::exhaustive::Search  & ptr) -> motif_search::SolutionOP {
                return ptr.next(); } )
            ;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// motif_search::monte_carlo
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // classes

    py::class_<motif_search::monte_carlo::DefaultScorer, std::shared_ptr<motif_search::monte_carlo::DefaultScorer>>(m, "MonteCarloDefaultScorer")
            // ctors
            .def(py::init<>())
                    // methods
            .def("clone",[] (motif_search::monte_carlo::DefaultScorer const & ptr) -> motif_search::monte_carlo::Scorer * {
                return ptr.clone(); } )
            .def("score",[] (motif_search::monte_carlo::DefaultScorer  & ptr, structure::BasepairState const & bps) ->  float {
                return ptr.score(bps); } )
                    // inherited methods
            .def("clone",[] (motif_search::monte_carlo::Scorer const & ptr) -> motif_search::monte_carlo::Scorer * {
                return ptr.clone(); } )
            .def("set_target",[] (motif_search::monte_carlo::Scorer  & ptr, structure::BasepairStateOP target, bool target_an_aligned_end) {
                ptr.set_target(target, target_an_aligned_end); } )
            .def("score",[] (motif_search::monte_carlo::Scorer  & ptr, structure::BasepairState const & bps) -> float {
                return ptr.score(bps); } )
            ;

    py::class_<motif_search::monte_carlo::MotifSwapMove, std::shared_ptr<motif_search::monte_carlo::MotifSwapMove>>(m, "MotifSwapMove")
            // ctors
            .def(py::init<motif_search::monte_carlo::ScorerOP, motif_search::SolutionToplogy const &>())
                    // methods
            .def("apply",[] (motif_search::monte_carlo::MotifSwapMove  & ptr, motif_data_structure::MotifStateGraphOP msg, float current_score) -> bool {
                return ptr.apply(msg, current_score); } )
            .def("score",[] (motif_search::monte_carlo::MotifSwapMove  & ptr) -> float {
                return ptr.score(); } )
            .def("undo",[] (motif_search::monte_carlo::MotifSwapMove  & ptr, motif_data_structure::MotifStateGraphOP msg) {
                ptr.undo(msg); } )
                    // inherited methods
            .def("apply",[] (motif_search::monte_carlo::Move  & ptr, motif_data_structure::MotifStateGraphOP msgop, float num) -> bool {
                return ptr.apply(msgop, num); } )
            .def("score",[] (motif_search::monte_carlo::Move  & ptr) -> float {
                return ptr.score(); } )
            .def("undo",[] (motif_search::monte_carlo::Move  & ptr, motif_data_structure::MotifStateGraphOP msgop) {
                ptr.undo(msgop); } )
            .def("set_temperature",[] (motif_search::monte_carlo::Move  & ptr, float temp) {
                ptr.set_temperature(temp); } )
            .def("scale_temperature",[] (motif_search::monte_carlo::Move  & ptr, float scale) {
                ptr.scale_temperature(scale); } )
            ;
/*
        py::class_<Move, std::shared_ptr<Move>>(m, "Move")
		// ctors
		.def(py::init<String const &>())
		// methods
		.def("apply",[] (Move  & ptr, motif_data_structure:: MotifStateGraphOP, float ) -> bool {
		 return ptr.apply(MotifStateGraphOP, ); } )
		.def("score",[] (Move  & ptr) -> float {
		 return ptr.score(); } )
		.def("undo",[] (Move  & ptr, motif_data_structure:: MotifStateGraphOP) {
		ptr.undo(MotifStateGraphOP); } )
		.def("set_temperature",[] (Move  & ptr, float temp) {
		ptr.set_temperature(temp); } )
		.def("scale_temperature",[] (Move  & ptr, float scale) {
		ptr.scale_temperature(scale); } )
		;
		;

        py::class_<Parameters, std::shared_ptr<Parameters>>(m, "Parameters")
		// public attributes
		.def_readwrite("accept_score", &Parameters::accept_score)
		.def_readwrite("max_size", &Parameters::max_size)
		;
*/
    py::class_< motif_search::monte_carlo::MoveSet, std::shared_ptr< motif_search::monte_carlo::MoveSet>>(m, "MoveSet");

    py::class_< motif_search::monte_carlo::ScaledScorer, std::shared_ptr< motif_search::monte_carlo::ScaledScorer>>(m, "ScaledScorer")
            // ctors
            .def(py::init<float,float>())
                    // methods
            .def("clone",[] ( motif_search::monte_carlo::ScaledScorer const & ptr) ->  motif_search::monte_carlo::Scorer * {
                return ptr.clone(); } )
            .def("score",[] ( motif_search::monte_carlo::ScaledScorer  & ptr, structure::BasepairState const & bps) ->  float {
                return ptr.score(bps); } )
                    // inherited methods
            .def("clone",[] ( motif_search::monte_carlo::Scorer const & ptr) ->  motif_search::monte_carlo::Scorer * {
                return ptr.clone(); } )
            .def("set_target",[] ( motif_search::monte_carlo::Scorer  & ptr, structure::BasepairStateOP target, bool target_an_aligned_end) {
                ptr.set_target(target, target_an_aligned_end); } )
            .def("score",[] ( motif_search::monte_carlo::Scorer  & ptr, structure::BasepairState const & bps) -> float {
                return ptr.score(bps); } )
            ;
/*
        py::class_<Scorer, std::shared_ptr<Scorer>>(m, "Scorer")
		// ctors
		.def(py::init<>())
		// methods
		.def("clone",[] (Scorer const & ptr) -> Scorer * {
		 return ptr.clone(); } )
		.def("set_target",[] (Scorer  & ptr, motif_search::monte_carlostructure::BasepairStateOP target, bool target_an_aligned_end) {
		ptr.set_target(target, target_an_aligned_end); } )
		.def("score",[] (Scorer  & ptr, structure::BasepairState const & bps) -> float {
		 return ptr.score(bps); } )
		;
*/
    py::class_<motif_search::monte_carlo::ScorerFactory, std::shared_ptr<motif_search::monte_carlo::ScorerFactory>>(m, "MonteCarloScorerFactory")
            // ctors
            .def(py::init<>())
                    // methods
            .def("get_scorer",[] (motif_search::monte_carlo::ScorerFactory  & ptr, String const & scorer_name) -> motif_search::monte_carlo::ScorerOP {
                return ptr.get_scorer(scorer_name); } )
            ;

    py::class_<motif_search::monte_carlo::Search, std::shared_ptr<motif_search::monte_carlo::Search>>(m, "MonteCarloSearch")
            // ctors
            .def(py::init<motif_search::monte_carlo::ScorerOP, motif_search::SolutionToplogy const &,motif_search::SolutionFilterOP>())
                    // methods
            .def("clone",[] (motif_search::monte_carlo::Search const & ptr) -> motif_search::Search * {
                return ptr.clone(); } )
            .def("setup",[] (motif_search::monte_carlo::Search  & ptr, motif_search::ProblemOP p) {
                ptr.setup(p); } )
            .def("start",[] (motif_search::monte_carlo::Search  & ptr) {
                ptr.start(); } )
            .def("finished",[] (motif_search::monte_carlo::Search  & ptr) -> bool {
                return ptr.finished(); } )
            .def("next",[] (motif_search::monte_carlo::Search  & ptr) -> motif_search::SolutionOP {
                return ptr.next(); } )
            ;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// motif_search::path_finding
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    m.def("default_selector", [] () -> SelectorOP {
        return default_selector(); }
    );
    // classes

    py::class_<AstarScorer, std::shared_ptr<AstarScorer>>(m, "AstarScorer")
            // ctors
            .def(py::init<>())
                    // methods
            .def("clone",[] (AstarScorer const & ptr) -> motif_search::path_finding::Scorer * {
                return ptr.clone(); } )
            .def("score",[] (AstarScorer  & ptr, motif_search::path_finding::Node const & node) ->  float {
                return ptr.score(node); } )
            .def("score",[] (AstarScorer  & ptr, motif::MotifState & ms, motif_search::path_finding::Node const & node) ->  float {
                return ptr.score(ms, node); } )
                    // inherited methods
            .def("clone",[] (motif_search::path_finding::Scorer const & ptr) -> motif_search::path_finding::Scorer * {
                return ptr.clone(); } )
            .def("set_target",[] (motif_search::path_finding::Scorer  & ptr, structure::BasepairStateOP target, bool target_an_aligned_end) {
                ptr.set_target(target, target_an_aligned_end); } )
            .def("score",[] (motif_search::path_finding::Scorer  & ptr, motif_search::path_finding::Node const & node) ->  float {
                return ptr.score(node); } )
            .def("score",[] (motif_search::path_finding::Scorer  & ptr, motif::MotifState & state, motif_search::path_finding::Node const & node) ->  float {
                return ptr.score(state, node); } )
            .def("accept_score",[] (motif_search::path_finding::Scorer  & ptr, motif_search::path_finding::Node const & node) ->  float {
                return ptr.accept_score(node); } )
            .def("set_dummy",[] (motif_search::path_finding::Scorer  & ptr, float dummy) {
                ptr.set_dummy(dummy); } )
            ;

    py::class_<GreedyScorer, std::shared_ptr<GreedyScorer>>(m, "GreedyScorer")
            // ctors
            .def(py::init<>())
                    // methods
            .def("clone",[] (GreedyScorer const & ptr) -> motif_search::path_finding::Scorer * {
                return ptr.clone(); } )
            .def("score",[] (GreedyScorer  & ptr, motif_search::path_finding::Node const & node) ->  float {
                return ptr.score(node); } )
            .def("score",[] (GreedyScorer  & ptr, motif::MotifState & ms, motif_search::path_finding::Node const & node) ->  float {
                return ptr.score(ms, node); } )
                    // inherited methods
            .def("clone",[] (motif_search::path_finding::Scorer const & ptr) -> motif_search::path_finding::Scorer * {
                return ptr.clone(); } )
            .def("set_target",[] (motif_search::path_finding::Scorer  & ptr, structure::BasepairStateOP target, bool target_an_aligned_end) {
                ptr.set_target(target, target_an_aligned_end); } )
            .def("score",[] (motif_search::path_finding::Scorer  & ptr, motif_search::path_finding::Node const & node) ->  float {
                return ptr.score(node); } )
            .def("score",[] (motif_search::path_finding::Scorer  & ptr, motif::MotifState & state, motif_search::path_finding::Node const & node) ->  float {
                return ptr.score(state, node); } )
            .def("accept_score",[] (motif_search::path_finding::Scorer  & ptr, motif_search::path_finding::Node const & node) ->  float {
                return ptr.accept_score(node); } )
            .def("set_dummy",[] (motif_search::path_finding::Scorer  & ptr, float dummy) {
                ptr.set_dummy(dummy); } )
            ;

    py::class_<motif_search::path_finding::Node, std::shared_ptr<motif_search::path_finding::Node>>(m, "MotifSearchPathFindingNode")
            // ctors
            .def(py::init<motif::MotifStateOP,motif_search::path_finding::NodeOP,float,int,int,int>())
                    // methods
            .def("level",[] (motif_search::path_finding::Node const & ptr) -> int {
                return ptr.level(); } )
            .def("size",[] (motif_search::path_finding::Node const & ptr) ->  int {
                return ptr.size(); } )
            .def("node_type",[] (motif_search::path_finding::Node const & ptr) ->  int {
                return ptr.node_type(); } )
            .def("ss_score",[] (motif_search::path_finding::Node const & ptr) ->  float {
                return ptr.ss_score(); } )
            .def("score",[] (motif_search::path_finding::Node const & ptr) ->  float {
                return ptr.score(); } )
            .def("parent_end_index",[] (motif_search::path_finding::Node const & ptr) ->  int {
                return ptr.parent_end_index(); } )
            .def("state",[] (motif_search::path_finding::Node const & ptr) ->  motif::MotifStateOP {
                return ptr.state(); } )
            .def("parent",[] (motif_search::path_finding::Node const & ptr) ->  motif_search::path_finding::NodeOP {
                return ptr.parent(); } )
            ;
/*
        py::class_<NodeCompare, std::shared_ptr<NodeCompare>>(m, "NodeCompare")
		// operators
		.def(py::self () motif_search::path_finding::NodeOP)
		;
*/


    py::class_<RoundRobinSelector, std::shared_ptr<RoundRobinSelector>>(m, "RoundRobinSelector")
            // ctors
            .def(py::init<>())
                    // methods
            .def("clone",[] (RoundRobinSelector const & ptr) -> Selector * {
                return ptr.clone(); } )
            .def("add",[] (RoundRobinSelector  & ptr, String const & lib_name) {
                ptr.add(lib_name); } )
            .def("add",[] (RoundRobinSelector  & ptr, String const & lib_name, motif::MotifStateEnsembleOP mse) {
                ptr.add(lib_name, mse); } )
            .def("add",[] (RoundRobinSelector  & ptr, motif::MotifOP motif) {
                ptr.add(motif); } )
                    // inherited methods
            .def("clone",[] (Selector const & ptr) -> Selector * {
                return ptr.clone(); } )
            .def("add",[] (Selector  & ptr, String const & lib_name) {
                ptr.add(lib_name); } )
            .def("add",[] (Selector  & ptr, String const & lib_name, motif::MotifStateEnsembleOP mse) {
                ptr.add(lib_name, mse); } )
            .def("add",[] (Selector  & ptr, motif::MotifOP motif) {
                ptr.add(motif); } )
            .def("connect",[] (Selector  & ptr, String const & name_i, String const & name_j) {
                ptr.connect(name_i, name_j); } )
            .def("start",[] (Selector const & ptr, int parent_type) {
                ptr.start(parent_type); } )
            .def("next",[] (Selector const & ptr) -> SelectorNodeDataOP {
                return ptr.next(); } )
            .def("finished",[] (Selector const & ptr) -> bool {
                return ptr.finished(); } )
            .def("size",[] (Selector const & ptr) -> size_t {
                return ptr.size(); } )
            ;
/*
        py::class_<Scorer, std::shared_ptr<Scorer>>(m, "Scorer")
		// ctors
		.def(py::init<>())
		// methods
		.def("clone",[] (Scorer const & ptr) -> Scorer * {
		 return ptr.clone(); } )
		.def("set_target",[] (Scorer  & ptr, motif_search::path_findingstructure::BasepairStateOP target, bool target_an_aligned_end) {
		ptr.set_target(target, target_an_aligned_end); } )
		.def("score",[] (Scorer  & ptr, Node const &) ->  float {
		 return ptr.score(&); } )
		.def("score",[] (Scorer  & ptr, motif::MotifState &, Node const &) ->  float {
		 return ptr.score(&, &); } )
		.def("accept_score",[] (Scorer  & ptr, Node const & node) ->  float {
		 return ptr.accept_score(node); } )
		.def("set_dummy",[] (Scorer  & ptr, float dummy) {
		ptr.set_dummy(dummy); } )
		;
*/
    py::class_<motif_search::path_finding::Search, std::shared_ptr<motif_search::path_finding::Search>>(m, "PathFindingSearch")
            // ctors
            .def(py::init<motif_search::path_finding::ScorerOP,SelectorOP,motif_search::SolutionFilterOP>())
                    // methods
            .def("clone",[] (motif_search::path_finding::Search const & ptr) -> motif_search::Search * {
                return ptr.clone(); } )
            .def("setup",[] (motif_search::path_finding::Search  & ptr, motif_search:: ProblemOP prblm) {
                ptr.setup(prblm); } )
            .def("start",[] (motif_search::path_finding::Search  & ptr) {
                ptr.start(); } )
            .def("finished",[] (motif_search::path_finding::Search  & ptr) -> bool {
                return ptr.finished(); } )
            .def("next",[] (motif_search::path_finding::Search  & ptr) -> motif_search::SolutionOP {
                return ptr.next(); } )
            ;

    py::class_<Selector, std::shared_ptr<Selector>>(m, "Selector")
            // ctors
            .def(py::init<>())
                    // methods
            .def("clone",[] (Selector const & ptr) -> Selector * {
                return ptr.clone(); } )
            .def("add",[] (Selector  & ptr, String const & lib_name) {
                ptr.add(lib_name); } )
            .def("add",[] (Selector  & ptr, String const & lib_name, motif::MotifStateEnsembleOP mse) {
                ptr.add(lib_name, mse); } )
            .def("add",[] (Selector  & ptr, motif::MotifOP motif) {
                ptr.add(motif); } )
            .def("connect",[] (Selector  & ptr, String const & name_i, String const & name_j) {
                ptr.connect(name_i, name_j); } )
            .def("start",[] (Selector const & ptr, int parent_type) {
                ptr.start(parent_type); } )
            .def("next",[] (Selector const & ptr) -> SelectorNodeDataOP {
                return ptr.next(); } )
            .def("finished",[] (Selector const & ptr) -> bool {
                return ptr.finished(); } )
            .def("size",[] (Selector const & ptr) -> size_t {
                return ptr.size(); } )
            ;

    py::class_<SelectorNodeData, std::shared_ptr<SelectorNodeData>>(m, "SelectorNodeData")
            // ctors
            .def(py::init<String const &,motif::MotifStateOPs const &,int>())
                    // public attributes
            .def_readwrite("name", &SelectorNodeData::name)
            .def_readwrite("motif_states", &SelectorNodeData::motif_states)
            .def_readwrite("type", &SelectorNodeData::type)
            ;

}