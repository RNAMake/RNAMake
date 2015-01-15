import ensemble_utils
import eterna_utils
import inv_utils
import sys
import re
import random
import os
import math
import simplejson
import subprocess


strategy_names = ['merryskies_only_as_in_the_loops', 'aldo_repetition', 'dejerpha_basic_test', 'eli_blue_line', 'clollin_gs_in_place', 'quasispecies_test_by_region_boundaries', 'eli_gc_pairs_in_junction', 'eli_no_blue_nucleotides_in_hook', 'mat747_31_loops', 'merryskies_1_1_loop', 'xmbrst_clear_plot_stack_caps_and_safe_gc', 'jerryp70_jp_stratmark', 'eli_energy_limit_in_tetraloops', 'eli_double_AUPair_strategy', 'eli_green_blue_strong_middle_half', 'eli_loop_pattern_for_small_multiloops', 'eli_tetraloop_similarity', 'example_gc60', 'penguian_clean_dotplot', 'eli_twisted_basepairs', 'aldo_loops_and_stacks', 'eli_direction_of_gc_pairs_in_multiloops_neckarea', 'eli_multiloop_similarity', 'eli_green_line', 'ding_quad_energy', 'quasispecies_test_by_region_loops', 'berex_berex_loop_basic', 'eli_legal_placement_of_GUpairs', 'merryskies_1_1_loop_energy', 'ding_tetraloop_pattern', 'aldo_mismatch', 'eli_tetraloop_blues', 'eli_red_line', 'eli_wrong_direction_of_gc_pairs_in_multiloops', 'deivad_deivad_strategy', 'eli_direction_of_gc_pairs_in_multiloops', 'eli_no_blue_nucleotides_strategy', 'berex_basic_test', 'eli_numbers_of_yellow_nucleotides_pr_length_of_string', 'kkohli_test_by_kkohli']
ensemble = ensemble_utils.Ensemble("L2", strategy_names, None)


designs = eterna_utils.get_synthesized_designs_from_eterna_server(False, "http://eterna.cmu.edu/eterna_get_synthesized.php?all=1")
csv = open("data.csv", "w+")

first = True
for design in designs:
    scoremap = ensemble.score(design)
    if first:
        header_str = "design name"
        for key in scoremap:
            if key != "finalscore" and key != "finalscore_normalized":
                header_str += ", " + key.replace(",", " ")
        header_str += ", chemical reactivity score"
        first= False
        csv.write("%s\n" % header_str)
    
    entry = design['soltitle'].replace(",", " ")
    for key in scoremap:
        if key != "finalscore" and key != "finalscore_normalized":
            if scoremap[key] > -90000:
                entry += ", %f" % scoremap[key]
            else:
                entry += ", N/A"
    
    if design['score'] > -90000:
        entry += ", %f" % design['score']
    else:
        entry += ", N/A"
    csv.write("%s\n" % entry)
    
csv.close()
