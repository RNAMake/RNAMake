import ensemble_utils
import ensemble_design
import settings

class SequenceDesigner:

    def __init__(self, op="sparse"):
        self._setup_ensemble(op)

    def _setup_ensemble(self, op):
        weights=None
        strategy_names = ['merryskies_only_as_in_the_loops', 'aldo_repetition',
            'dejerpha_basic_test', 'eli_blue_line', 'clollin_gs_in_place',
            'quasispecies_test_by_region_boundaries', 'eli_gc_pairs_in_junction',
            'eli_no_blue_nucleotides_in_hook', 'mat747_31_loops', 'merryskies_1_1_loop',
            'xmbrst_clear_plot_stack_caps_and_safe_gc', 'jerryp70_jp_stratmark',
            'eli_energy_limit_in_tetraloops', 'eli_double_AUPair_strategy',
            'eli_green_blue_strong_middle_half','eli_loop_pattern_for_small_multiloops',
            'eli_tetraloop_similarity', 'example_gc60', 'penguian_clean_dotplot',
            'eli_twisted_basepairs', 'aldo_loops_and_stacks',
            'eli_direction_of_gc_pairs_in_multiloops_neckarea','eli_multiloop_similarity',
            'eli_green_line', 'ding_quad_energy','quasispecies_test_by_region_loops',
            'berex_berex_loop_basic','eli_legal_placement_of_GUpairs',
            'merryskies_1_1_loop_energy', 'ding_tetraloop_pattern', 'aldo_mismatch',
            'eli_tetraloop_blues', 'eli_red_line',
            'eli_wrong_direction_of_gc_pairs_in_multiloops', 'deivad_deivad_strategy',
            'eli_direction_of_gc_pairs_in_multiloops', 'eli_no_blue_nucleotides_strategy',
            'berex_basic_test', 'eli_numbers_of_yellow_nucleotides_pr_length_of_string',
            'kkohli_test_by_kkohli']

        if op == "sparse":
            weights_file_name = "no_validation_training/weights_sparse_5.overall.txt"
            weights_f = open(settings.base_dir + "/" + weights_file_name,"r")
            weights = []
            for line in weights_f:
                weights.append(float(line))

        #for i, name in enumerate(strategy_names):
            #if weights[i] != 0:
            #    print name, weights[i]

        self.ensemble = ensemble_utils.Ensemble("sparse", strategy_names, weights)

    def design(self, structure, sequence, count=1):
        score_cutoff = 90
        if   len(structure) <= 50:
            score_cutoff = 70
        elif len(structure) <= 80:
            score_cutoff = 80

        solutions = []
        for i in range(count):
            res = ensemble_design.inverse_fold_whole(structure, sequence,
                                                     self.ensemble.score,
                                                     score_cutoff, "sparse")
            solutions.append(res)
        return solutions

if __name__ == "__main__":
    designer = SequenceDesigner()
    solutions = designer.design("(((....)))", "GGGUUCGCCC")
    #solutions = designer.design("(((((....)))))", "GGAGCUUCGGCUCC")
    #solutions = designer.design("(((....)))", "NNNUUCGNNN")
    print solutions

