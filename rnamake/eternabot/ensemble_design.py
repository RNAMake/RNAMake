from . import ensemble_utils
from . import eterna_utils
from . import inv_utils
import sys
import re
import random
import os
import math
import simplejson
import subprocess


def initial_sequence_with_gc_caps(secstruct, sequence_constraints, no_gccap):
    """
    args:
    secstruct is input secondary structure
    sequence_constraints are the constraints on the sequence
    no_gccap is boolean indicating whether to include gc cap

    return:
    string containing the initial sequence
    """
    n = len(secstruct)
    pairmap = eterna_utils.get_pairmap_from_secstruct(secstruct)
    sequence = []

    # initialize to all A except constrained positions
    for ii in range(0, n):
        if(sequence_constraints[ii] == "N"):
            sequence.append("A")
        else:
            sequence.append(sequence_constraints[ii])

    # put GC caps at end of pairs
    pair_count = 0
    gc_count = 0
    for ii in range(0, n):

        # skip constrained nucleotides
        if(sequence_constraints[ii] != "N"):
            if (pairmap[ii] > ii):
                pair_count += 1
            continue

        # identify caps
        if(pairmap[ii] > ii):
            pair_count += 1
            is_cap = False

            if(ii == 0):
                is_cap = True

            # identify nucleotides next to pairs
            if(ii > 0 and pairmap[ii-1] < 0):
                is_cap = True
            if(ii < n-1 and pairmap[ii+1] < 0):
                is_cap = True

            # JEEFIX for conventional
            #is_cap = False

            # randomly assign a pair
            if(is_cap and not no_gccap):
                if(random.random() > 0.5):
                    sequence[ii] = "G"
                    sequence[pairmap[ii]] = "C"
                else:
                    sequence[ii] = "C"
                    sequence[pairmap[ii]] = "G"

                gc_count += 1

    # modify the number of GC pairs to 0.6 of total
    max_gc = pair_count * 0.6
    gc_prob = float(max_gc - gc_count) / pair_count

    for ii in range(0, n):
        if(sequence_constraints[ii] != "N"):
            continue

        if(pairmap[ii] > ii):
            if(sequence[ii] != "G" and sequence[ii] != "C"):
                if(random.random() < gc_prob):
                    if(random.random() > 0.5):
                        sequence[ii] = "G"
                        sequence[pairmap[ii]] = "C"
                    else:
                        sequence[ii] = "C"
                        sequence[pairmap[ii]] = "G"
                    gc_count += 1
                else:
                    if(random.random() > 0.5):
                        sequence[ii] = "A"
                        sequence[pairmap[ii]] = "U"
                    else:
                        sequence[ii] = "U"
                        sequence[pairmap[ii]] = "A"

    # convert the sequence to a string
    string_sequence = ""
    for ii in range(0, n):
        string_sequence += sequence[ii]

    return string_sequence


def get_random_base():
    """
    generates random base

    args:
    none

    returns:
    random base as char
    """
    #bases = "GAUGCG"
    bases = "GCAU"
    randn = int(math.floor(random.random() * 4) % 4)
    return bases[randn]


def get_random_pair():
    """
    generates random base pair

    args:
    none

    return:
    random base pair as string
    """
    bases = ["GC", "CG", "AU", "UA"]
    randn = int(math.floor(random.random() * 4) % 4)
    return bases[randn]


def get_sequence_array(sequence):
    """
    turns sequence string into sequence array

    args:
    sequence is a string

    returns:
    sequence in array form
    """
    arr = []
    for ii in range(0, len(sequence)):
        arr.append(sequence[ii])

    return arr


def get_sequence_string(arr):
    """
    turns sequence array into string

    args:
    sequence is an array

    returns:
    sequence in array form
    """
    string = ""
    for ii in range(0, len(arr)):
        string += arr[ii]

    return string


def inverse_fold(
        secstruct,
        start_sequence,
        constraints,
        scoring_func,
        score_cutoff):
    """
    args:
    secstruct contains secondary structure
    start_sequence is the current sequence
    constraints contain constraints on sequence positions
    scoring_func is a function for scoring sequences
    score_cutoff is the cutoff score to stop the search

    return:
    best sequence, bp distance, final design
    """

    # if secstruct is all unpaired, then keep start sequence
    dotonly = True
    for ii in range(0, len(secstruct)):
        if(secstruct[ii] != "."):
            dotonly = False

    if(dotonly or len(secstruct) < 20):
        return [start_sequence, 0, 0, {}]

    bases = "GAUC"
    pairs = ["GC", "CG", "AU", "UA"]

    # print "%s %s" % (secstruct, start_sequence)

    target_pairmap = eterna_utils.get_pairmap_from_secstruct(secstruct)
    n = len(start_sequence)

    # get indices for stacks and loops
    stack_indices = []
    loop_indices = []
    for ii in range(0, len(target_pairmap)):
        if target_pairmap[ii] > ii and constraints[ii] == "N":
            stack_indices.append(ii)
        elif target_pairmap[ii] < 0 and constraints[ii] == "N":
            if(ii > 0 and target_pairmap[ii-1] >= 0):
                loop_indices.append(ii)
            elif(ii < len(target_pairmap)-1 and target_pairmap[ii+1] >= 0):
                loop_indices.append(ii)

    # get information about current sequence
    sequence = start_sequence
    native = inv_utils.fold(sequence)[0]
    bp_distance = eterna_utils.bp_distance(secstruct, native)
    native_pairmap = eterna_utils.get_pairmap_from_secstruct(native)
    design = eterna_utils.get_design_from_sequence(sequence, secstruct)
    design_score = scoring_func(design)

    # initialize variables for iteration
    walk_iter = 0
    stale_move = 0

    best_sequence = sequence
    best_bp_distance = bp_distance
    best_native = native
    best_native_pairmap = native_pairmap
    best_design = design
    best_design_score = design_score

    # get indices for positions that can change
    index_array = []
    for ii in range(0, n):
        if(constraints[ii] == "N"):
            index_array.append(ii)

    # loop as long as bp distance too large or design score too small
    while(bp_distance > 10 or design_score['finalscore'] < score_cutoff) and len(index_array) > 0:
        # print "%d %d %f" %(walk_iter, best_bp_distance, best_design_score['finalscore'])
        # random.shuffle(index_array)

        # iterate over nucleotides in sequence
        moved = False
        for rrr in range(0, len(index_array)):
            ii = index_array[rrr]

            break_for = False
            if(target_pairmap[ii] != native_pairmap[ii]):
                # nonmatching nucleotides that should be unpaired
                if(target_pairmap[ii] < 0):
                    # try each other possible base in this position
                    for jj in range(0, len(bases)):
                        if sequence[ii] == bases[jj]:
                            continue
                        mut_array = get_sequence_array(sequence)
                        mut_array[ii] = bases[jj]

                        mut_sequence = get_sequence_string(mut_array)
                        mut_native = inv_utils.fold(mut_sequence)[0]
                        mut_bp_distance = eterna_utils.bp_distance(
                            secstruct,
                            mut_native)
                        mut_design = eterna_utils.get_design_from_sequence(
                            mut_sequence,
                            secstruct)
                        mut_score = scoring_func(design)

                        # if distance or score is better for mutant, update the
                        # current sequence
                        if((mut_bp_distance < bp_distance) or mut_score['finalscore'] > design_score['finalscore']):
                            sequence = mut_sequence
                            native = mut_native
                            bp_distance = mut_bp_distance
                            native_pairmap = eterna_utils.get_pairmap_from_secstruct(
                                native)
                            design = mut_design
                            design_score = mut_score

                            # if distance or score is better for mutant than
                            # best, update the best sequence
                            if((mut_bp_distance < best_bp_distance or mut_bp_distance < 10 or len(secstruct) < 50) and mut_score['finalscore'] > best_design_score['finalscore']):
                                best_sequence = mut_sequence
                                best_bp_distance = mut_bp_distance
                                best_native = mut_native
                                best_native_pairmap = eterna_utils.get_pairmap_from_secstruct(
                                    mut_native)
                                best_design = mut_design
                                best_design_score = mut_score

                            break_for = True
                            break
                # nonmatching nucleotides that should be paired
                else:
                    # try each other possible pair in this position
                    current_pair = sequence[ii] + sequence[target_pairmap[ii]]
                    for jj in range(0, len(pairs)):
                        if current_pair == pairs[jj]:
                            continue
                        mut_array = get_sequence_array(sequence)
                        mut_array[ii] = pairs[jj][0]
                        mut_array[target_pairmap[ii]] = pairs[jj][1]

                        mut_sequence = get_sequence_string(mut_array)
                        mut_native = inv_utils.fold(mut_sequence)[0]
                        mut_bp_distance = eterna_utils.bp_distance(
                            secstruct,
                            mut_native)
                        mut_design = eterna_utils.get_design_from_sequence(
                            mut_sequence,
                            secstruct)
                        mut_score = scoring_func(design)

                        # if distance or score is better for mutant, update the
                        # current sequence
                        if((mut_bp_distance < bp_distance) or mut_score['finalscore'] > design_score['finalscore']):
                            sequence = mut_sequence
                            bp_distance = mut_bp_distance
                            native = mut_native
                            native_pairmap = eterna_utils.get_pairmap_from_secstruct(
                                native)
                            design = mut_design
                            design_score = mut_score

                            # if distance or score is better for mutant than
                            # best, update best
                            if((mut_bp_distance < best_bp_distance or mut_bp_distance < 10 or len(secstruct) < 50) and mut_score['finalscore'] > best_design_score['finalscore']):
                                best_sequence = mut_sequence
                                best_bp_distance = mut_bp_distance
                                best_native = mut_native
                                best_native_pairmap = eterna_utils.get_pairmap_from_secstruct(
                                    mut_native)
                                best_design = mut_design
                                best_design_score = mut_score

                            break_for = True
                            break
            # if the sequence has updated, break loop through nucleotides
            if(break_for):
                moved = True
                break

        # if sequence hasn't been updated, randomize sequence
        if(moved == False):

            stale_move += 1

            if(stale_move > 10):
                return [
                    best_sequence,
                    best_bp_distance,
                    best_design_score['finalscore'],
                    best_design_score]

            rand = random.random()

            mut_array = get_sequence_array(sequence)
            if(rand < 0.5 and len(stack_indices) > 0):
                rindex = stack_indices[
                    int(random.random() * len(stack_indices)) % len(stack_indices)]
                if(target_pairmap[rindex] < 0):
                    print "Something is wrong"
                    sys.exit(0)

                if(rand < 0.33):
                    temp = mut_array[rindex]
                    mut_array[rindex] = mut_array[target_pairmap[rindex]]
                    mut_array[target_pairmap[rindex]] = temp
                else:
                    pair = get_random_pair()
                    mut_array[rindex] = pair[0]
                    mut_array[target_pairmap[rindex]] = pair[1]
            else:
                if (len(loop_indices) > 0):
                    rindex = loop_indices[
                        int(random.random() * len(loop_indices)) % len(loop_indices)]
                    mut_array[rindex] = get_random_base()

            sequence = get_sequence_string(mut_array)
            native = inv_utils.fold(sequence)[0]
            bp_distance = eterna_utils.bp_distance(secstruct, native)
            design = eterna_utils.get_design_from_sequence(sequence, secstruct)
            design_score = scoring_func(design)
        else:
            stale_move = 0

        moved == False
        walk_iter += 1

        # if it has been too many iterations, finish
        if(walk_iter > 400):
            return [
                best_sequence,
                best_bp_distance,
                best_design_score['finalscore'],
                best_design_score]

    return [
        best_sequence,
        best_bp_distance,
        best_design_score['finalscore'],
        best_design_score]


def inverse_fold_whole(
        secstruct,
        constraints,
        score_func,
        score_cutoff,
        conventional):
    """
    args:
    secstruct contains secondary structure
    constraints contain constraints on sequence positions
    score_func is a function for scoring sequences, takes in dict of designs
        returns scores
    score_cutoff is the cutoff score to stop the search
    conventional is the method

    return:
    res is dictionary of results
    """
    # initialize sequence and corresponding folds/bp distances/designs
    start_sequence = initial_sequence_with_gc_caps(
        secstruct,
        constraints,
        conventional)

    start_native = inv_utils.fold(start_sequence)[0]
    start_bp = eterna_utils.bp_distance(start_native, secstruct)
    start_design = eterna_utils.get_design_from_sequence(
        start_sequence,
        secstruct)

    # store initial data to result dict
    res = {}
    res['initial'] = [start_sequence, start_bp, score_func(start_design)]

    # determine if N is present in constraints
    flag = 0

    for i in range(0, len(constraints)):
        if constraints[i] == 'N':
            flag = 1
            break

    # if not, only one possible sequence.  otherwise, recursively fold
    if flag == 0:
        res['end'] = [start_sequence, start_bp, score_func(start_design)]
    else:
        start_sequence_array = get_sequence_array(start_sequence)
        elements = eterna_utils.get_rna_elements_from_secstruct(secstruct)
        root = elements[0]
        inverse_fold_recursive(
            root,
            secstruct,
            start_sequence_array,
            constraints,
            score_func,
            score_cutoff)

        end_sequence = get_sequence_string(start_sequence_array)
        end_design = eterna_utils.get_design_from_sequence(
            end_sequence,
            secstruct)
        end_native = inv_utils.fold(end_sequence)[0]
        end_bp = eterna_utils.bp_distance(end_native, secstruct)

        res['end'] = [end_sequence, end_bp, score_func(end_design)]
    return res


def inverse_fold_recursive(
        root,
        secstruct,
        sequence_array,
        constraints,
        score_func,
        score_cutoff):
    """
    args:
    root is the root element
    secstruct contains secondary structure
    sequence_array contains the current sequence as an array
    constraints contain constraints on sequence positions
    score_func is a function for scoring sequences
    score_cutoff is the cutoff score to stop the search

    return:
    none
    """
    if(len(root.indices_) == 0):
        return

    if(len(root.children_) == 0):
        return

    # run this recursively over all children of each element
    for ii in range(0, len(root.children_)):
        inverse_fold_recursive(
            root.children_[ii],
            secstruct,
            sequence_array,
            constraints,
            score_func,
            score_cutoff)

    # do this for each element
    sequence = get_sequence_string(sequence_array)
    root_start_index = root.indices_[0]
    if(root.type_ == eterna_utils.RNAELEMENT_LOOP):
        root_end_index = root.indices_[len(root.indices_)-1]
    else:
        root_end_index = root.indices_[1]

    # print "RR %d %d %s" % (root_start_index, root_end_index,
    # str(root.indices_))

    res = inverse_fold(
        secstruct[
            root_start_index:root_end_index +
            1],
        sequence[
            root_start_index:root_end_index +
            1],
        constraints[
            root_start_index:root_end_index +
            1],
        score_func,
        score_cutoff)
    res_sequence = res[0]

    for ii in range(root_start_index, root_end_index+1):
        sequence_array[ii] = res_sequence[ii-root_start_index]


def main():
    # JEEFIX for conventional
    strategy_names = [
        'merryskies_only_as_in_the_loops', 'aldo_repetition',
        'dejerpha_basic_test', 'eli_blue_line', 'clollin_gs_in_place',
        'quasispecies_test_by_region_boundaries', 'eli_gc_pairs_in_junction',
        'eli_no_blue_nucleotides_in_hook', 'mat747_31_loops',
        'merryskies_1_1_loop', 'xmbrst_clear_plot_stack_caps_and_safe_gc',
        'jerryp70_jp_stratmark', 'eli_energy_limit_in_tetraloops',
        'eli_double_AUPair_strategy', 'eli_green_blue_strong_middle_half',
        'eli_loop_pattern_for_small_multiloops', 'eli_tetraloop_similarity',
        'example_gc60', 'penguian_clean_dotplot', 'eli_twisted_basepairs',
        'aldo_loops_and_stacks',
        'eli_direction_of_gc_pairs_in_multiloops_neckarea',
        'eli_multiloop_similarity', 'eli_green_line', 'ding_quad_energy',
        'quasispecies_test_by_region_loops', 'berex_berex_loop_basic',
        'eli_legal_placement_of_GUpairs', 'merryskies_1_1_loop_energy',
        'ding_tetraloop_pattern', 'aldo_mismatch', 'eli_tetraloop_blues',
        'eli_red_line', 'eli_wrong_direction_of_gc_pairs_in_multiloops',
        'deivad_deivad_strategy', 'eli_direction_of_gc_pairs_in_multiloops',
        'eli_no_blue_nucleotides_strategy', 'berex_basic_test',
        'eli_numbers_of_yellow_nucleotides_pr_length_of_string',
        'kkohli_test_by_kkohli']

    # parse inputs and specify parameters
    score_cutoff = 90

    secstructs = [sys.argv[3]]
    constraints = [sys.argv[4]]
    repeat = int(sys.argv[5])
    puzzle_titles = ["FMN Aptamer"]

    op = sys.argv[1]

    if len(secstructs) <= 50:
        score_cutoff = 70
    elif len(secstructs) <= 80:
        score_cutoff = 80

    if op == "L2":
        weights_file_name = "no_validation_training/weights_L2.overall.txt"
        scores_file_name = "no_validation_training/predicted_score_L2.overall.unnormalized.txt"
    elif op == "sparse":
        weights_file_name = "no_validation_training/weights_sparse_5.overall.txt"
        scores_file_name = "no_validation_training/predicted_score_sparse_5.overall.unnormalized.txt"
        weights_f = open(weights_file_name, "r")
        weights = []
        for line in weights_f:
            weights.append(float(line))

    elif op == "conventional":
        strategy_names = [
            'example_gc60',
            'penguian_clean_dotplot',
            'berex_simplified_berex_test']
        weights_file_name = "no_validation_training/weights_conventional.overall.txt"
        scores_file_name = "no_validation_training/predicted_score_L2.overall.unnormalized.txt"
    else:
        print "No type specified"
        sys.exit(0)

    # create ensemble of different strategies
    is_sparse_test = (op == "sparse")
    if(is_sparse_test):
        ensemble = ensemble_utils.Ensemble(op, strategy_names, weights)
    else:
        ensemble = ensemble_utils.Ensemble(op, strategy_names, None)

    # loop over all secondary structures
    for pp in range(0, len(secstructs)):
        outputs = []
        print "\n\n======%s==========\n\n" % puzzle_titles[pp]

        # repeat folding according to argument
        for ii in range(0, repeat):
            outputs.append(
                (inverse_fold_whole(
                    secstructs[pp],
                    constraints[pp],
                    ensemble.score,
                    score_cutoff,
                    op == "conventional")))

        # sort outputs based on final score of designs
        outputs = sorted(
            outputs,
            key=lambda output: output['end'][2]['finalscore'],
            reverse=True)
        outputs.append(secstructs[pp])

        init_sum = 0
        final_sum = 0

        # write output
        print "printing output"
        outputfile = open(
            "/persistent/drupal/html/files/" +
            sys.argv[2] +
            ".txt",
            "w+")
        outputfile.write(simplejson.dumps(outputs))
        outputfile.close()
        print "printed output"


if __name__ == "__main__":
    main()
