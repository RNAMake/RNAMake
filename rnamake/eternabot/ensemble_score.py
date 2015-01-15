import ensemble_utils
import eterna_utils
import sys
import re



#JEEFIX for conventional
strategy_names = ['merryskies_only_as_in_the_loops', 'aldo_repetition', 'dejerpha_basic_test', 'eli_blue_line', 'clollin_gs_in_place', 'quasispecies_test_by_region_boundaries', 'eli_gc_pairs_in_junction', 'eli_no_blue_nucleotides_in_hook', 'mat747_31_loops', 'merryskies_1_1_loop', 'xmbrst_clear_plot_stack_caps_and_safe_gc', 'jerryp70_jp_stratmark', 'eli_energy_limit_in_tetraloops', 'eli_double_AUPair_strategy', 'eli_green_blue_strong_middle_half', 'eli_loop_pattern_for_small_multiloops', 'eli_tetraloop_similarity', 'example_gc60', 'penguian_clean_dotplot', 'eli_twisted_basepairs', 'aldo_loops_and_stacks', 'eli_direction_of_gc_pairs_in_multiloops_neckarea', 'eli_multiloop_similarity', 'eli_green_line', 'ding_quad_energy', 'quasispecies_test_by_region_loops', 'berex_berex_loop_basic', 'eli_legal_placement_of_GUpairs', 'merryskies_1_1_loop_energy', 'ding_tetraloop_pattern', 'aldo_mismatch', 'eli_tetraloop_blues', 'eli_red_line', 'eli_wrong_direction_of_gc_pairs_in_multiloops', 'deivad_deivad_strategy', 'eli_direction_of_gc_pairs_in_multiloops', 'eli_no_blue_nucleotides_strategy', 'berex_basic_test', 'eli_numbers_of_yellow_nucleotides_pr_length_of_string', 'kkohli_test_by_kkohli']
score_cutoff = 90

secstructs =  [".....(((((((((((((((....)))))))((((......((((....)))).....))))(((((((....)))))))))))))))...................."]
constraints = ["GGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGGAUAUNNNNNNNNNNAGAAGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAAGAAACAACAACAACAAC"]
puzzle_titles = ["FMN Aptamer"]

op = sys.argv[1]

if op == "L2":
	weights_file_name = "no_validation_training/weights_L2.overall.txt"
	scores_file_name = "no_validation_training/predicted_score_L2.overall.unnormalized.txt"
elif op == "sparse":
	weights_file_name = "no_validation_training/weights_sparse_5.overall.txt"
	scores_file_name = "no_validation_training/predicted_score_sparse_5.overall.unnormalized.txt"
elif op == "conventional":
	strategy_names = ['example_gc60', 'penguian_clean_dotplot', 'berex_simplified_berex_test']
	weights_file_name = "no_validation_training/weights_conventional.overall.txt"
	scores_file_name = "no_validation_training/predicted_score_L2.overall.unnormalized.txt"
else:
	print "No type specified"
	sys.exit(0)

weight_file = open(weights_file_name,'r')
weights = []
iter = 0
for line in weight_file:
	if(iter != 26):
		weights.append(float(line))
	iter += 1
weight_file.close()

scores_file = open(scores_file_name,'r')
scores = []
iter = 0
for line in scores_file:
	if(iter > 0):
		scores.append(float(line))
	iter += 1
scores_file.close()

is_sparse_test = op == "sparse"
if(is_sparse_test):
	ensemble = ensemble_utils.Ensemble(op, strategy_names, weights)
else:
	ensemble = ensemble_utils.Ensemble(op, strategy_names, None)

ensemble.test_weights(weights)
ensemble.test_scores(scores)

all_designs = eterna_utils.get_synthesized_designs_from_eterna_server(False, "http://eterna.cmu.edu/eterna_get_synthesized.php?all=2")
for design in all_designs:
	try:
		puztitle = design['puztitle']
		soltitle = design['soltitle']
		score = ensemble.score(design)['finalscore']
		print "%s %s %f" % (puztitle,soltitle,score)
