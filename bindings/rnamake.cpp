// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/cast.h>
#include <pybind11/operators.h>

// base files
#include <base/application.hpp>
#include <base/global_constants.h>
#include <base/backtrace.h>
#include <base/settings.h>
#include <base/command_line_parser.hpp>
#include <base/file_io.h>
#include <base/env_manager.h>
#include <base/sys_interface.h>
#include <base/vector_container.h>
#include <base/string.h>
#include <base/assertions.h>
#include <base/exception.h>
#include <base/util.hpp>
#include <base/log.h>
#include <base/cl_option.h>
#include <base/option.h>
#include <base/types.h>

// data_structure files
#include <data_structure/graph.h>
#include <data_structure/graph_adjacency_list.h>
#include <data_structure/graph_iter_list.h>
#include <data_structure/graph_base.h>

// data_structure/graph files
#include <data_structure/graph/graph.h>
#include <data_structure/graph/graph_node.h>

// data_structure/tree files
#include <data_structure/tree/tree_node.h>
#include <data_structure/tree/tree.h>

// eternabot files
#include <eternabot/feature_generator.h>
#include <eternabot/sequence_designer.h>
#include <eternabot/scorer.h>
#include <eternabot/strategy.h>

// eternabot/strategy files
#include <eternabot/strategy/num_of_yellow.h>
#include <eternabot/strategy/direction_of_gc.h>
#include <eternabot/strategy/clear_plot.h>
#include <eternabot/strategy/berex_test.h>
#include <eternabot/strategy/a_basic_test.h>

// main files

// math files
#include <math/xyz_vector.h>
#include <math/quaternion.h>
#include <math/xyz_matrix.h>
#include <math/hashing.h>
#include <math/numerical.h>
#include <math/transform.h>
#include <math/euler.h>
#include <math/stats.h>

// motif files
#include <motif/motif_state_aligner.h>
#include <motif/motif_state_ensemble.h>
#include <motif/motif.h>
#include <motif/motif_scorer.h>
#include <motif/motif_state.h>
#include <motif/motif_ensemble.h>
//#include <motif/pose_factory.h>
//#include <motif/pose.h>
#include <motif/motif_to_secondary_structure.h>
#include <motif/motif_factory.h>

// motif_data_structure files
#include <motif_data_structure/motif_graph.h>
#include <motif_data_structure/motif_state_tree.h>
#include <motif_data_structure/motif_tree.h>
#include <motif_data_structure/motif_state_ensemble_graph.h>
#include <motif_data_structure/motif_connection.h>
#include <motif_data_structure/motif_state_graph.hpp>
#include <motif_data_structure/motif_topology.h>
#include <motif_data_structure/motif_state_node.hpp>
#include <motif_data_structure/motif_state_ensemble_tree.h>
#include <motif_data_structure/motif_merger.h>

// motif_search files
#include <motif_search/problem.h>
#include <motif_search/solution_filter.h>
#include <motif_search/solution_topology.h>
#include <motif_search/search.h>
#include <motif_search/motif_state_monte_carlo.h>

// motif_search/exhaustive files
#include <motif_search/exhaustive/scorer.h>
#include <motif_search/exhaustive/search.h>
#include <motif_search/exhaustive/motif_state_enumerator.h>

// motif_search/monte_carlo files
#include <motif_search/monte_carlo/scorer.h>
#include <motif_search/monte_carlo/search.h>

// motif_search/path_finding files
#include <motif_search/path_finding/node.h>
#include <motif_search/path_finding/scorer.h>
#include <motif_search/path_finding/selector.h>
#include <motif_search/path_finding/search.h>

// motif_tools files
#include <motif_tools/segmenter.h>

// motif_tree_state_alt_pather files
//#include <motif_tree_state_alt_pather.h>

// resources files
#include <resources/motif_ensemble_sqlite_connection.h>
#include <resources/motif_state_sqlite_library.h>
#include <resources/motif_sqlite_library.h>
#include <resources/motif_sqlite_connection.h>
#include <resources/added_motif_library.h>
#include <resources/resource_manager.h>
#include <resources/sqlite_library.h>
#include <resources/motif_state_ensemble_sqlite_library.h>

// secondary_structure files
#include <secondary_structure/secondary_structure_tree.h>
#include <secondary_structure/motif.h>
#include <secondary_structure/pose.h>
#include <secondary_structure/sequence_tools.h>
#include <secondary_structure/secondary_structure_parser.h>
#include <secondary_structure/residue.h>
#include <secondary_structure/basepair.h>
#include <secondary_structure/util.h>
#include <secondary_structure/chain.h>
#include <secondary_structure/rna_structure.h>
#include <secondary_structure/sequence_constraint.h>
#include <secondary_structure/structure.h>

// sequence_optimization files
#include <sequence_optimization/sequence_optimizer.h>
#include <sequence_optimization/sequence_optimizer_3d.hpp>

// structure files
#include <structure/structure.h>
#include <structure/pdb_parser.h>
#include <structure/close_chain.h>
#include <structure/chain.h>
#include <structure/exceptions.h>
#include <structure/basepair_state.h>
#include <structure/beads.h>
#include <structure/residue_type_set_manager.h>
#include <structure/residue_type.h>
#include <structure/residue_type_set.h>
#include <structure/cif_parser.h>
#include <structure/atom.h>
#include <structure/residue.h>
#include <structure/basepair.h>
#include <structure/is_equal.h>
#include <structure/rna_structure.h>

// thermo_fluctuation files
#include <thermo_fluctuation/thermo_fluc_simulation.h>
#include <thermo_fluctuation/thermo_fluc_relax.h>
//#include <thermo_fluctuation/thermo_fluc_simulation_devel.h>
#include <thermo_fluctuation/thermo_fluc_sampler.h>
#include <thermo_fluctuation/thermo_fluc_scorer.h>

// thermo_fluctuation/graph files
#include <thermo_fluctuation/graph/simulation.h>
#include <thermo_fluctuation/graph/sterics.h>
#include <thermo_fluctuation/graph/scorer.h>
#include <thermo_fluctuation/graph/logging.h>
#include <thermo_fluctuation/graph/sampler.h>

// util files
#include <util/motif_type.h>
#include <util/monte_carlo.h>
#include <util/steric_lookup.hpp>
#include <util/cartesian_product.h>
#include <util/x3dna.h>
#include <util/basic_io.hpp>
#include <util/uuid.h>
#include <util/dssr.h>
#include <util/sqlite3_connection.h>
#include <util/random_number_generator.h>
#include <util/csv.h>

// vienna files
#include <vienna/energy_par.h>
#include <vienna/vienna.h>
#include <vienna/pair_mat.h>


// Bindings start

namespace py = pybind11;

PYBIND11_MODULE(pyRNAMake, m) {

// Free functions


/*
    m.def("motif_data_structure_motif_state_ensemble_graph_from_motif_graph", &motif_data_structure::motif_state_ensemble_graph_from_motif_graph,
	py::arg("&"),
	py::arg("&"),
	py::arg("&"),
	py::arg("&"));


    m.def("motif_data_structure_graph_to_tree", &motif_data_structure::graph_to_tree,
	py::arg("mg"),
	py::arg("start") = nullptr,
	py::arg("last_end") = nullptr);
        

    m.def("resources_get_motif_from_resource_manager", &resources::get_motif_from_resource_manager,
	py::arg("name") = resources::dummy_name,
	py::arg("end_id") = resources::dummy_end_id,
	py::arg("end_name") = resources::dummy_name);
        

    m.def("resources_build_sqlite_library", &resources::build_sqlite_library,
	py::arg("path"),
	py::arg("data"),
	py::arg("keys"),
	py::arg("primary_key"));


//m.def("resources_sqlite3_escape", &resources::sqlite3_escape,
//#py::arg("unescaped_string"));
        
    m.def("util_type_to_str", &util::type_to_str,
	py::arg("mtype"));
        

    m.def("util_str_to_type", &util::str_to_type,
	py::arg("s"));
        

    m.def("util_get_x3dna_by_type", &util::get_x3dna_by_type,
	py::arg("name"));
        

    m.def("util_get_str_from_x3dna_type", &util::get_str_from_x3dna_type,
	py::arg("type"));
        

    m.def("util_compare_bps", &util::compare_bps,
	py::arg("lhs"),
	py::arg("rhs"));
        

    m.def("util_json_cleanup", &util::json_cleanup);
        
    m.def("util_points_to_pdb_str", &util::points_to_pdb_str,
	py::arg("points"));
        

    m.def("util_points_to_pdb", &util::points_to_pdb,
	py::arg("filename"),
	py::arg("points"));
        

    m.def("util_get_double", &util::get_double,
	py::arg("&"),
	py::arg("&"));
        

    m.def("util_get_int", &util::get_int,
	py::arg("&"),
	py::arg("&"));
        

    m.def("util_get_char", &util::get_char,
	py::arg("&"),
	py::arg("&"));
        

    m.def("util_get_string", &util::get_string,
	py::arg("&"),
	py::arg("&"));
        

    m.def("util_get_point", &util::get_point,
	py::arg("&"),
	py::arg("&"));
        

    m.def("util_get_quaternion", &util::get_quaternion,
	py::arg("&"),
	py::arg("&"));
        

    m.def("util_get_matrix", &util::get_matrix,
	py::arg("&"));
        

    m.def("util_get_reals", &util::get_reals,
	py::arg("&"),
	py::arg("&"));
        

    m.def("util_get_ints", &util::get_ints,
	py::arg("&"),
	py::arg("&"));
        

    m.def("util_get_nts", &util::get_nts,
	py::arg("&"));
        

    m.def("util_get_pairs", &util::get_pairs,
	py::arg("&"));
        

    m.def("util_get_hairpins", &util::get_hairpins,
	py::arg("&"));
        

    m.def("util_get_helices", &util::get_helices,
	py::arg("&"));
        

    m.def("util_get_stems", &util::get_stems,
	py::arg("&"));
        

    m.def("util_get_iloops", &util::get_iloops,
	py::arg("&"));
        

    m.def("util_get_juncs", &util::get_juncs,
	py::arg("&"));
        

    m.def("util_get_elements", &util::get_elements,
	py::arg("&"),
	py::arg("&"),
	py::arg("&"),
	py::arg("&"),
	py::arg("&"),
	py::arg("&"),
	py::arg("&"));
        
    m.def("io_detail_chop_next_column", &io::detail::chop_next_column,
	py::arg("line"),
	py::arg("col_begin"),
	py::arg("col_end"));
        

    m.def("io_detail_parse_line", &io::detail::parse_line,
	py::arg("line"),
	py::arg("sorted_col"),
	py::arg("col_order"));
        

    m.def("io_detail_parse_header_line", &io::detail::parse_header_line,
	py::arg("line"),
	py::arg("col_order"),
	py::arg("col_name"),
	py::arg("ignore_policy"));
        

    m.def("io_detail_parse", &io::detail::parse,
	py::arg("col"),
	py::arg("x"));
        

    m.def("io_detail_parse_unsigned_integer", &io::detail::parse_unsigned_integer,
	py::arg("col"),
	py::arg("x"));
        

    m.def("io_detail_parse_signed_integer", &io::detail::parse_signed_integer,
	py::arg("col"),
	py::arg("x"));


    m.def("io_detail_parse_float", &io::detail::parse_float,
	py::arg("col"),
	py::arg("x"));

    m.def("motif_align_motif", &motif::align_motif,
	py::arg("ref_bp_state"),
	py::arg("motif_end"),
	py::arg("motif"));
        

    m.def("motif_get_aligned_motif", &motif::get_aligned_motif,
	py::arg("ref_bp"),
	py::arg("motif_end"),
	py::arg("motif"));
        

    m.def("motif_ref_motif", &motif::ref_motif);
        

    m.def("motif_file_to_motif", &motif::file_to_motif,
	py::arg("path"));
        

    m.def("motif_clash_between_motifs", &motif::clash_between_motifs,
	py::arg("m1"),
	py::arg("m2"),
	py::arg("clash_radius") = 2.7);
        

    m.def("motif_align_motif_state", &motif::align_motif_state,
	py::arg("ref_bp_state"),
	py::arg("org_state"));
        

    m.def("motif_get_aligned_motif_state", &motif::get_aligned_motif_state,
	py::arg("ref_bp_state"),
	py::arg("cur_state"),
	py::arg("org_state"));
        

    m.def("math_vector_from_str", &math::vector_from_str,
	py::arg("s"));
        

    m.def("math_vector_to_str", &math::vector_to_str,
	py::arg("v"));
        

    m.def("math_vectors_to_str", &math::vectors_to_str,
	py::arg("vs"));
        

    m.def("math_vectors_from_str", &math::vectors_from_str,
	py::arg("s"));
        

    m.def("math_get_random_quaternion", &math::get_random_quaternion);
        

    m.def("math_get_quaternion_from_matrix", &math::get_quaternion_from_matrix,
	py::arg("m"));
        

    m.def("math_power_iteration", &math::power_iteration,
	py::arg("A"),
	py::arg("eigen_values"),
	py::arg("num_simulations"));
        

    m.def("math_dot_vector", [] ( math::Matrix const& m, math::Vector const& v, math::Vector & vr )  { math::dot_vector(m, v, vr);},
	py::arg("m"),
	py::arg("v"),
	py::arg("vr"));
        

    m.def("math_dot_vectors", &math::dot_vectors,
	py::arg("m"),
	py::arg("v"),
	py::arg("vr"));
        

    m.def("math_matrix_from_str", &math::matrix_from_str,
	py::arg("s"));
        

    m.def("math_transform_1", [] (math::Matrix& m1) -> math::Matrix { return math::transform_1(m1);},
	py::arg("m"));
        

    m.def("math_matrix_to_str", &math::matrix_to_str,
	py::arg("m"));
        

    m.def("math_are_floats_equal", &math::are_floats_equal,
	py::arg("a"),
	py::arg("b"),
	py::arg("tol") = 0.001);
        

    m.def("math_are_xyzVector_equal", &math::are_xyzVector_equal,
	py::arg("vec"),
	py::arg("correct_vec"),
	py::arg("tol") = 0.001);
        

    m.def("math_are_xyzVectors_equal", &math::are_xyzVectors_equal,
	py::arg("v"),
	py::arg("vc"));
        

    m.def("math_are_xyzMatrix_equal", &math::are_xyzMatrix_equal,
	py::arg("m"),
	py::arg("mc"));
        

    //m.def("math_roughly_equal", &math::roughly_equal,
	//py::arg("v1"),
	//py::arg("v2"),
	//py::arg("tolerance") = 0.001);
        

    //m.def("math_>", &math::>,
	//py::arg("v1"),
	//py::arg("v2"),
	//py::arg("tolerance"));
        

    m.def("math_norm", &math::norm,
	py::arg("v"));
        

    m.def("math_calc_euler", &math::calc_euler,
	py::arg("M"),
	py::arg("euler"));
        

    m.def("math_axis_angle_from_matrix", &math::axis_angle_from_matrix,
	py::arg("m"),
	py::arg("aa"));
        

    m.def("math_degrees", &math::degrees,
	py::arg("radians"));
        

    m.def("math_sum", &math::sum,
	py::arg("a"));
        

    m.def("math_sqsum", &math::sqsum,
	py::arg("a"));
        

    m.def("math_stdev", &math::stdev,
	py::arg("nums"));
        

    m.def("math_mean", &math::mean,
	py::arg("a"));
        

    m.def("math_pearson_coeff", &math::pearson_coeff,
	py::arg("x"),
	py::arg("y"));
        

    m.def("math_avg_unsigned_diff", &math::avg_unsigned_diff,
	py::arg("x"),
	py::arg("y"));
        

    m.def("structure_create_coord_system", &structure::create_coord_system,
	py::arg("atoms"));
        

    m.def("structure_to_radians", &structure::to_radians,
	py::arg("degrees"));
        

    m.def("structure_virtual_atom", &structure::virtual_atom,
	py::arg("name"),
	py::arg("l"),
	py::arg("theta"),
	py::arg("phi"),
	py::arg("parent_atoms"));
        

    m.def("structure_get_projection", &structure::get_projection,
	py::arg("coord"),
	py::arg("current_pos"),
	py::arg("projection_axis"));
        

    m.def("structure_axis_angle_to_rot_matrix", &structure::axis_angle_to_rot_matrix,
	py::arg("angle"),
	py::arg("al"));
        

    m.def("structure_close_torsion", &structure::close_torsion,
	py::arg("which_dir"),
	py::arg("parent_atoms"),
	py::arg("daughter_atoms"),
	py::arg("match_atoms_1"),
	py::arg("match_atoms_2"));
        

    m.def("structure_get_res_ref_frame", &structure::get_res_ref_frame,
	py::arg("r"));
        

    m.def("structure_replace_missing_phosphate_backbone", &structure::replace_missing_phosphate_backbone,
	py::arg("r"),
	py::arg("r_template"));
        

    m.def("structure_close_chain", &structure::close_chain,
	py::arg("chain"));
        

    m.def("structure_str_to_basepairstate", &structure::str_to_basepairstate,
	py::arg("s"));
        

    m.def("structure_get_ref_bp_state", &structure::get_ref_bp_state);
        

    m.def("structure_get_bpstate_rotation_diff", &structure::get_bpstate_rotation_diff,
	py::arg("bp1"),
	py::arg("bp2"));
        

    m.def("structure_frame_distance", &structure::frame_distance,
	py::arg("current"),
	py::arg("end"),
	py::arg("endflip"));
        

    m.def("structure_new_score_function", &structure::new_score_function,
	py::arg("current"),
	py::arg("end"),
	py::arg("endflip"));
        

    m.def("structure_new_score_function_new", &structure::new_score_function_new,
	py::arg("current"),
	py::arg("end"),
	py::arg("endflip"));
        

    m.def("structure_are_basepair_states_equal", &structure::are_basepair_states_equal,
	py::arg("a"),
	py::arg("b"));
        

    m.def("structure_center", &structure::center,
	py::arg("atoms"));
        

    m.def("structure_wc_bp", &structure::wc_bp,
	py::arg("bp"));
        

    m.def("structure_gu_bp", &structure::gu_bp,
	py::arg("bp"));
        

    m.def("structure_are_atoms_equal", &structure::are_atoms_equal,
	py::arg("a1"),
	py::arg("a2"),
	py::arg("tol") = 0.001);
        

    m.def("structure_are_atom_vectors_equal", &structure::are_atom_vectors_equal,
	py::arg("atoms_1"),
	py::arg("atoms_2"),
	py::arg("tol") = 0.001);
        

    m.def("structure_are_residues_equal", &structure::are_residues_equal,
	py::arg("r1"),
	py::arg("r2"),
	py::arg("check_uuids") = 1);
        

    m.def("structure_are_chains_equal", &structure::are_chains_equal,
	py::arg("c1"),
	py::arg("c2"),
	py::arg("check_uuids") = 1);
        

    m.def("structure_are_structures_equal", &structure::are_structures_equal,
	py::arg("s1"),
	py::arg("s2"),
	py::arg("check_uuids") = 1);
        
    m.def("structure_m_dot_v", &structure::m_dot_v,
	py::arg("m"),
	py::arg("v"));
        

    m.def("structure_connect_residues_into_chains", &structure::connect_residues_into_chains,
	py::arg("residues"),
	py::arg("chains"));
        

    m.def("structure_end_from_basepairs", &structure::end_from_basepairs,
	py::arg("s"),
	py::arg("bps"));
        

    m.def("structure_subselect_basepairs_with_res", &structure::subselect_basepairs_with_res,
	py::arg("res"),
	py::arg("all_bps"));

//    m.def("secondary_structure_tree_from_pose", &secondary_structure::tree_from_pose,
//	py::arg("p"));
        

    m.def("secondary_structure_get_res_types_from_sequence", &secondary_structure::get_res_types_from_sequence,
	py::arg("sequence"),
	py::arg("residue_types"));
        

//    m.def("secondary_structure_find_res_types_in_pose", &secondary_structure::find_res_types_in_pose,
//	py::arg("p"),
//	py::arg("residue_types"));
        

    m.def("secondary_structure_find_gc_helix_stretches", &secondary_structure::find_gc_helix_stretches,
	py::arg("p"),
	py::arg("length"));
        

    m.def("secondary_structure_find_longest_gc_helix_stretch", &secondary_structure::find_longest_gc_helix_stretch,
	py::arg("p"));
        

    m.def("secondary_structure_convert_res_name_to_type", &secondary_structure::convert_res_name_to_type,
	py::arg("c"));
        

    m.def("secondary_structure_is_gc_pair", &secondary_structure::is_gc_pair,
	py::arg("bp"));
        

    m.def("secondary_structure_is_au_pair", &secondary_structure::is_au_pair,
	py::arg("bp"));
        

    m.def("secondary_structure_is_gu_pair", &secondary_structure::is_gu_pair,
	py::arg("bp"));
        

    m.def("secondary_structure_assign_end_id", &secondary_structure::assign_end_id,
	py::arg("ss"),
	py::arg("end"));
        

    m.def("secondary_structure_fill_basepairs_in_ss", &secondary_structure::fill_basepairs_in_ss,
	py::arg("ss"));
        

    m.def("motif_search_path_finding_default_selector", &motif_search::path_finding::default_selector);
        

    m.def("base_get_os_name", &base::get_os_name);
        

    m.def("base_resources_path", &base::resources_path);
        

    m.def("base_lib_path", &base::lib_path);
        

    m.def("base_motif_dirs", &base::motif_dirs);
        

    m.def("base_x3dna_path", &base::x3dna_path);
        

    m.def("base_unittest_resource_dir", &base::unittest_resource_dir);
        

    m.def("base_demangle", &base::demangle,
	py::arg("string"));
        

    m.def("base_print_backtrace", &base::print_backtrace);
        

    m.def("base_save_backtrace", &base::save_backtrace);
        

    m.def("base_file_exists", &base::file_exists,
	py::arg("name"));
        

    m.def("base_is_dir", &base::is_dir,
	py::arg("path"));
        

    m.def("base_get_lines_from_file", &base::get_lines_from_file,
	py::arg("noexcept(false"));
        
    
    m.def("base_execute_command", &base::execute_command,
	py::arg("cmd"));
        

    m.def("base_execute_command_json", &base::execute_command_json,
	py::arg("cmd"));

    m.def("base_split_str_by_delimiter", &base::split_str_by_delimiter,
	py::arg("s"),
	py::arg("delimiter"));
        

    m.def("base_tokenize_line", &base::tokenize_line,
	py::arg("raw_line"));
        

    m.def("base_join_by_delimiter", &base::join_by_delimiter,
	py::arg("strs"),
	py::arg("delimiter"));
        

    m.def("base_filename", &base::filename,
	py::arg("path"));
        

    //m.def("base_base_dir", &base::base_dir,
	//py::arg("path"));
        

    m.def("base_determine_string_data_type", &base::determine_string_data_type,
	py::arg("str"));
        

    m.def("base_is_number", &base::is_number,
	py::arg("s"));
        

    m.def("base_ltrim", &base::ltrim,
	py::arg("s"));
        

    m.def("base_rtrim", &base::rtrim,
	py::arg("s"));
        

    m.def("base_trim", &base::trim,
	py::arg("s"));
        

    m.def("base_replace_all", &base::replace_all,
	py::arg("context"),
	py::arg("from"),
	py::arg("to"));
        

    m.def("expects", &expects,
      //py::arg("condition"),
      //py::arg("message"));
    

    m.def("ensures", &ensures,
      //py::arg("condition"),
      //py::arg("message"));
    

    m.def("element_in_vector", &element_in_vector,
	//py::arg("element"),
	//py::arg("vec"));
        

    m.def("base_init_logging", &base::init_logging,
	py::arg("log_level"));
        

    m.def("base_log_level_from_str", &base::log_level_from_str,
	py::arg("s"));
*/

}
