
#include "design_rna_scaffold/design_rna_scaffold.h"

#ifndef DESIGNRNASCAFFOLD_OPTIONS_H
#define DESIGNRNASCAFFOLD_OPTIONS_H

String
valid_pdb (String &path) {
    auto ending = path.substr(path.size() - 4);
    return ending == ".pdb" ? String{""} : String{"the file specified by --pdb must end in .pdb"};
}

String
valid_bp (String &bp) {
    const auto bp_pattern = std::regex("\b[A-Z][0-9]*-[A-Z][0-9]*\b");
    auto sm = std::smatch{};
    std::regex_match(bp, sm, bp_pattern);

    return sm.size() == 1 ? String{""} : String{bp + " is an invalid bp format"};
}

struct Parameters {

	struct Core{
		String pdb = "";
		String start_bp = "";
		String end_bp = "";
		int designs = 1;
		String log_level = "info";
		String extra_pdbs = "";
	};

	struct IO{
		bool dump_intermediate_pdbs = "";
		bool dump_pdbs = false;
		bool dump_scaffold_pdbs = false;
		String new_ensembles_file = "";
		bool no_out_file = false;
		String out_file = "default.out";
		String score_file = "default.scores";
	};

	struct Search{
		String ending_helix = "";
		String exhaustive_scorer = "default";
		int max_helix_length = 99;
		String mc_scorer = "default";
		int min_helix_length = 4;
		String motif_path = "";
		bool no_basepair_checks = false;
		float scaled_score_d = 1.0f;
		float scaled_score_r = 2.0f;
		float cutoff = 7.5f;
		int max_size = 999999;
		String type = "path_finding";
		String solution_filter = "RemoveDuplicateHelices";
		String starting_helix = "";
		bool no_sterics = false;
		bool only_tether_opt = false;
		int max_motifs = 999;
	};

	struct SequenceOpt{
		bool skip = false;
		int sequences_per_design = 1;
		int steps = 10000;
	};

	struct ThermoFluc{
		bool perform = false;
		int steps = 1000000;
	};


	Core core = Core();
	IO io = IO();
	Search search = Search();
	SequenceOpt seq_opt = SequenceOpt();
	ThermoFluc thermo_fluc = ThermoFluc();

};

    Parameters parameters_ = Parameters();

void
DesignRNAScaffold::setup_options () {

	app_.add_option_group("core");
	app_.add_option_group("io");
	app_.add_option_group("search");
	app_.add_option_group("seq_opt");
	app_.add_option_group("thermo_fluc");
	app_.add_option("--pdb", parameters_.core.pdb, "path to a PDB file with input RNA structure")
		->default_val("")
		->group("core");

	app_.add_option("--start_bp", parameters_.core.start_bp, "starting basepair to be used in structure format: [CHAIN ID][NT1 NUM]-[CHAIN ID][NT2 NUM]")
		->default_val("")
		->group("core");

	app_.add_option("--end_bp", parameters_.core.end_bp, "ending basepair to be used in structure format: [CHAIN ID][NT1 NUM]-[CHAIN ID][NT2 NUM]")
		->default_val("")
		->group("core");

	app_.add_option("--designs", parameters_.core.designs, "number of designs to create. Default is 1")
		->default_val(1)
		->group("core");

	app_.add_option("--log_level", parameters_.core.log_level, "level for global logging")
		->default_val("info")
		->group("core");

	app_.add_option("--extra_pdbs", parameters_.core.extra_pdbs, "deliminted list of other pdbs used in building")
		->default_val("")
		->group("core");

	app_.add_option("--dump_intermediate_pdbs", parameters_.io.dump_intermediate_pdbs, "flag to dump intermediate pdbs")
		->default_val("")
		->group("io");

	app_.add_flag("--dump_pdbs", parameters_.io.dump_pdbs, "TODO")
		->default_val(false)
		->group("io");

	app_.add_flag("--dump_scaffold_pdbs", parameters_.io.dump_scaffold_pdbs, "flag to output pdbs of just the design scaffold WITHOUT initial RNA structure useful for really big structures")
		->default_val(false)
		->group("io");

	app_.add_option("--new_ensembles", parameters_.io.new_ensembles_file, "flag to include new structural ensembles")
		->default_val("")
		->group("io");

	app_.add_flag("--no_out_file", parameters_.io.no_out_file, "if you only want the summary and not the actual structures")
		->default_val(false)
		->group("io");

	app_.add_option("--out_file", parameters_.io.out_file, "output file that contains all information to rebuild solutions")
		->default_val("default.out")
		->group("io");

	app_.add_option("--score_file", parameters_.io.score_file, "name of output file containining scoring information for design")
		->default_val("default.scores")
		->group("io");

	app_.add_option("--ending_helix", parameters_.search.ending_helix, "ending helix for design solution. Format = [TODO]")
		->default_val("")
		->group("search");

	app_.add_option("--exhaustive_scorer", parameters_.search.exhaustive_scorer, "TODO")
		->default_val("default")
		->group("search");

	app_.add_option("--max_helix_length", parameters_.search.max_helix_length, "maximum number of basepairs in a solution helix")
		->default_val(99)
		->group("search");

	app_.add_option("--mc_scorer", parameters_.search.mc_scorer, "TODO")
		->default_val("default")
		->group("search");

	app_.add_option("--min_helix_length", parameters_.search.min_helix_length, "minimum number of basepairs in a solution helix")
		->default_val(4)
		->group("search");

	app_.add_option("--motif_path", parameters_.search.motif_path, "TBD")
		->default_val("")
		->group("search");

	app_.add_flag("--no_basepair_checks", parameters_.search.no_basepair_checks, "flag to disable basepair checks")
		->default_val(false)
		->group("search");

	app_.add_option("--scaled_score_d", parameters_.search.scaled_score_d, "TODO")
		->default_val(1.0f)
		->group("search");

	app_.add_option("--scaled_score_r", parameters_.search.scaled_score_r, "TODO")
		->default_val(2.0f)
		->group("search");

	app_.add_option("--search_cutoff", parameters_.search.cutoff, "TODO")
		->default_val(7.5f)
		->group("search");

	app_.add_option("--search_max_size", parameters_.search.max_size, "maximum number of steps for a design search")
		->default_val(999999)
		->group("search");

	app_.add_option("--search_type", parameters_.search.type, "search type for traversing motif space")
		->default_val("path_finding")
		->group("search");

	app_.add_option("--solution_filter", parameters_.search.solution_filter, "TODO")
		->default_val("RemoveDuplicateHelices")
		->group("search");

	app_.add_option("--starting_helix", parameters_.search.starting_helix, "starting helix for design solution. Format = [TODO]")
		->default_val("")
		->group("search");

	app_.add_flag("--no_sterics", parameters_.search.no_sterics, "turns off sterics checks againsts supplied RNA structure")
		->default_val(false)
		->group("search");

	app_.add_flag("--only_tether_opt", parameters_.search.only_tether_opt, "gnore supplied structure other than sterics")
		->default_val(false)
		->group("search");

	app_.add_option("--search_max_motifs", parameters_.search.max_motifs, "TODO")
		->default_val(999)
		->group("search");

	app_.add_flag("--skip_sequence_optimization", parameters_.seq_opt.skip, "flag to skip sequence optimization of the design")
		->default_val(false)
		->group("seq_opt");

	app_.add_option("--sequences_per_design", parameters_.seq_opt.sequences_per_design, "number of sequences to try per motif design")
		->default_val(1)
		->group("seq_opt");

	app_.add_option("--seq_opt_steps", parameters_.seq_opt.steps, "TODO")
		->default_val(10000)
		->group("seq_opt");

	app_.add_flag("--thermo_fluc", parameters_.thermo_fluc.perform, "run thermo fluc procedure to estimate thermo fluc of helices")
		->default_val(false)
		->group("thermo_fluc");

	app_.add_option("--thermo_fluc_steps", parameters_.thermo_fluc.steps, "TODO")
		->default_val(1000000)
		->group("thermo_fluc");


}

#endif //DESIGNRNASCAFFOLD_OPTIONS_H
