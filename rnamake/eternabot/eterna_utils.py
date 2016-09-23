import csv
import sys
import random
import math
import vienna_parameters
import os
import re
import imp
import settings

import rnamake.settings


UNSCORABLE = -99999

RNAELEMENT_LOOP = "LOOP"
RNAELEMENT_STACK = "STACK"
DEFAULT_TEMPERATURE = 37.0

class RNAElement:
    def __init__(self):
        self.indices_ = []
        self.type_ = -1
        self.branching_stacks_ = 0
        self.score_ = 0
        self.parent_ = False
        self.children_ = []
        self.quad_scores_ = []

    def get_loop_groups(self):
        if(self.type_ != RNAELEMENT_LOOP):
            return []

        last_index = -999
        last_group = []
        groups = []
        for ii in range(0,len(self.indices_)):
            if(self.indices_[ii] != last_index+1 and last_index >= 0):
                groups.append(last_group)
                last_group = []

            last_group.append(self.indices_[ii])
            last_index = self.indices_[ii]

        if(len(last_group) > 0):
            groups.append(last_group)

        return groups

    def get_stack_length(self):
        if(self.type_ != RNAELEMENT_STACK):
            return 0
        return len(self.indices_)/2

    def get_pair_from_stack(self,pair_index,sequence):
        pair = sequence[self.indices_[pair_index * 2]] + sequence[self.indices_[pair_index * 2 + 1]]
        return pair.upper()

    def get_loop_closing_pairs(self,sequence,pairmap):
        if(self.type_ != RNAELEMENT_LOOP):
            return []

        pairs = []
        pair_indices = []

        if(self.parent_):
            parent = self.parent_
            if(parent.type_ != RNAELEMENT_STACK):
                print "ERROR Loop's parent is not a stack"
                sys.exit(0)

            npi = len(parent.indices_)
            pair = sequence[parent.indices_[npi-2]] + sequence[parent.indices_[npi-1]]
            pairs.append(pair)

        for ii in range(0,len(self.children_)):
            child = self.children_[ii]
            if(child.type_ != RNAELEMENT_STACK):
                print "ERROR Loop's child is not a stack"
                sys.exit(0)

            nci = len(child.indices_)
            pair = sequence[child.indices_[0]] + sequence[child.indices_[1]]
            pairs.append(pair)

        return pairs


def get_feature_vector_from_secstruct(secstruct):
    elements = get_rna_elements_from_secstruct(secstruct)
    pairmap = get_pairmap_from_secstruct(secstruct)

    dummy_sequence = ""
    for ii in range(0,len(secstruct)):
        dummy_sequence += "A"

    num_stacks = 0
    min_length_stack = -1
    max_length_stack = -1
    min_length_loop = -1
    max_length_loop = -1
    num_hairpins = 0
    num_multiloops = 0
    num_internal_loops = 0
    num_bulges = 0


    for ii in range(0,len(elements)):
        if(elements[ii].type_ == RNAELEMENT_STACK):
            num_stacks += 1
            slen = len(elements[ii].indices_)/2
            if(slen <min_length_stack or min_length_stack < 0):
                min_length_stack = slen
            if(slen >max_length_stack or max_length_stack < 0):
                max_length_stack = slen
        else:
            loop_groups = elements[ii].get_loop_groups()
            closing_pairs = elements[ii].get_loop_closing_pairs(dummy_sequence,pairmap)

            num_lg = len(loop_groups)
            num_cp = len(closing_pairs)

            if(not (0 in elements[ii].indices_)):
                llen = len(elements[ii].indices_)
                if(llen < min_length_loop or min_length_loop < 0):
                    min_length_loop = llen
                if(llen > max_length_loop or max_length_loop < 0):
                    max_length_loop = llen

            if(num_lg == 1 and num_cp == 1):
                num_hairpins += 1
            elif(num_lg == 1 and num_cp == 2):
                num_bulges += 1
            elif(num_lg == 2 and num_cp == 1):
                num_multiloops += 1
            elif(num_lg == 2 and num_cp == 2):
                num_internal_loops += 1
            elif(num_lg > 2):
                num_multiloops += 1
            elif(num_cp > 2):
                num_multiloops += 1
            else:
                print "ERROR : unknown loop type with %d loop groups and %d closing pairs" % (num_lg, num_cp)
                print "ERROR : " + str(loop_groups)

    return [num_stacks, min_length_stack, max_length_stack, min_length_loop, max_length_loop, num_hairpins, num_multiloops, num_internal_loops, num_bulges]

def get_pairmap_from_secstruct(secstruct):
    """
    generates dictionary containing pair mappings

    args:
    secstruct contains secondary structure string

    returns:
    dictionary with pair mappings
    """
    pair_stack = []
    pairs_array = []
    i_range = range(0,len(secstruct))

    # initialize all values to -1, meaning no pair
    for ii in i_range:
        pairs_array.append(-1)

    # assign pairs based on secstruct
    for ii in i_range:
        if(secstruct[ii] == "("):
            pair_stack.append(ii)
        elif(secstruct[ii] == ")"):
            index = pair_stack.pop()
            pairs_array[index] = ii
            pairs_array[ii] = index

    return pairs_array

def bp_distance(secstruct1, secstruct2):
    """
    calculates bp distance between two secondary structures

    args:
    secstruct1 is the first secondary structure
    secstruct2 is the second secondary structure

    returns:
    bp distance between structures
    """
    # ensure that secondary structures are the same length
    if(len(secstruct1) != len(secstruct2)):
        print "SS1 and SS2 lengths don't match"
        sys.exit(0)

    # generate pair mappings
    pairmap1 = get_pairmap_from_secstruct(secstruct1)
    pairmap2 = get_pairmap_from_secstruct(secstruct2)

    # +1 for each pair or single that doesn't match
    dist = 0

    for ii in range(0,len(pairmap1)):
        if(pairmap1[ii] != pairmap2[ii]):
            if(pairmap1[ii] > ii):
                dist += 1
            if(pairmap2[ii] > ii):
                dist += 1

    return dist

def get_rna_elements_from_secstruct(secstruct):
    """
    get array of pairs and stacks

    args:
    secstruct is the secondary structure

    return:
    list of RNA elements
    """
    elements = []
    pair_stack = []
    pairs_array = []
    i_range = range(0,len(secstruct))

    for ii in i_range:
        pairs_array.append(-1)

    for ii in i_range:
        if(secstruct[ii] == "("):
            pair_stack.append(ii)
        elif(secstruct[ii] == ")"):
            index = pair_stack.pop()
            pairs_array[index] = ii
            pairs_array[ii] = index

    get_rna_elements_from_secstruct_recursive(pairs_array,0,len(secstruct)-1,elements,-999,-999, False)
    return elements

def get_rna_elements_from_secstruct_recursive(pairs_array,start_index,end_index,elements, last_pair_start, last_pair_end, last_parent):
    """
    args:
    pairs_array is array containing the index of the pair of each nucleotide, negative i    f no pair
    start_index represents the begining of the current sequence
    end_index represents the end of the current sequence
    elements is a list of RNA elements
    last_pair_start contains the start index of the last pair
    last_pair_end contains the end index of the last pair
    last_parent is most recent parent element
    """

    # create a loop element
    new_element = RNAElement()
    new_element.type_ = RNAELEMENT_LOOP

    ii = start_index
    while( ii <= end_index):
        # unpaired nucleotides
        if(pairs_array[ii] < 0):
            new_element.indices_.append(ii)

            # create new stack element including all nucleotides since last pair
            if(last_pair_start >= 0):
                stack_element = RNAElement()
                stack_element.type_ = RNAELEMENT_STACK
                stack_element.parent_ = last_parent
                if(last_parent):
                    last_parent.children_.append(stack_element)

                for jj in range(last_pair_start,ii):
                    stack_element.indices_.append(jj)
                    stack_element.indices_.append(pairs_array[jj])

                elements.append(stack_element)
                last_pair_start = -999
                last_pair_end = -999
                new_element.parent_ = stack_element
                stack_element.children_.append(new_element)

            ii += 1
        # paired nucleotides
        elif(ii < pairs_array[ii]):
            # no unpaired since start of this iteration
            if(ii == start_index and pairs_array[ii] == end_index):
                if(last_pair_start < 0):
                    last_pair_start = ii
                    last_pair_end = pairs_array[ii]
                    last_parent = new_element
                get_rna_elements_from_secstruct_recursive(pairs_array,ii+1,pairs_array[ii]-1,elements,last_pair_start,last_pair_end, last_parent)
                return
            # at least one unpaired since start of this iteration
            else:
                # create stack element including all nucleotides since last pair
                if(last_pair_start >= 0):
                    stack_element = RNAElement()
                    stack_element.type_ = RNAELEMENT_STACK
                    stack_element.parent_ = last_parent
                    if(last_parent):
                        last_parent.children_.append(stack_element)

                    for jj in range(last_pair_start,ii):
                        stack_element.indices_.append(jj)
                        stack_element.indices_.append(pairs_array[jj])

                    elements.append(stack_element)
                    new_element.parent_ = stack_element
                    stack_element.children_.append(new_element)
                    last_pair_start = -999
                    last_pair_end = -999
                # recursively run for every pair
                get_rna_elements_from_secstruct_recursive(pairs_array,ii+1,pairs_array[ii]-1,elements,ii,pairs_array[ii], new_element)
            # skip parts of sequence covered by recursion
            ii = pairs_array[ii]+1
        else:
            print(str(start_index) + " " + str(end_index))
            print("ERROR : should not get here " + str(ii) + " " + str(pairs_array[ii]))
            sys.exit(0)

    #if(len(new_element.indices_) > 0):
    elements.append(new_element)

def validate_sequence_secstruct(sequence,secstruct):
    if(len(sequence) != len(secstruct)):
        return False

    pairmap = get_pairmap_from_secstruct(secstruct)

    for ii in range(0,len(pairmap)):
        if(pairmap[ii] > ii):
            pair = sequence[ii] + sequence[pairmap[ii]]

            if(pair != "GC" and pair != "CG" and pair != "AU" and pair != "UA" and pair != "GU" and pair !="UG"):
                return False

    return True

def get_GC_ratio(sequence,secstruct):
    pairmap = get_pairmap_from_secstruct(secstruct)

    total_pairs = 0
    gc_pairs = 0

    for ii in range(0,len(pairmap)):
        if(pairmap[ii] > ii):
            pair = sequence[ii] + sequence[pairmap[ii]]

            total_pairs += 1

            if(pair == "GC" or pair == "CG"):
                gc_pairs += 1

    return float(gc_pairs)/total_pairs



def fill_energy(elements,sequence,pairmap):
    for ii in range(0,len(elements)):
        if(elements[ii].type_ == RNAELEMENT_LOOP):
            loop_groups = elements[ii].get_loop_groups();
            num_loop_groups = len(loop_groups)
            closing_pairs = elements[ii].get_loop_closing_pairs(sequence,pairmap)
            num_closing_pairs = len(closing_pairs)

            if(elements[ii].indices_.count(0) > 0 or elements[ii].indices_.count(len(sequence) -1) >0):
                elements[ii].score_ = vienna_parameters.ml_energy(pairmap,sequence,0,True) / 100.0
            else:
                if(num_closing_pairs == 1 and num_loop_groups == 1):
                    hp_start = loop_groups[0][0]
                    hp_end = loop_groups[0][len(loop_groups[0])-1]
                    pair_type = vienna_parameters.pair_type(sequence[hp_start-1],sequence[hp_end+1])
                    size = hp_end - hp_start + 1
                    elements[ii].score_ = vienna_parameters.hairpin_energy(size, pair_type, vienna_parameters.letter_to_sequence_type(sequence[hp_start]), vienna_parameters.letter_to_sequence_type(sequence[hp_end]), sequence, hp_start-1,hp_end+1) / 100.0

                elif(num_closing_pairs == 2):
                    if(num_loop_groups > 2):
                        print "ERROR 2 closing pairs but 3 or more loop groups"
                        sys.exit(0)

                    n1 = 0
                    n2 = 0

                    if(num_loop_groups == 1):
                        n1 = len(loop_groups[0])
                        pi = loop_groups[0][0]-1
                        pj = pairmap[pi]

                        if(pi < pj):
                            i = loop_groups[0][0]-1
                            j = pairmap[i]

                            p = loop_groups[0][n1-1]+1
                            q = pairmap[p]
                        else:
                            q = loop_groups[0][0]-1
                            p = pairmap[q]

                            j = loop_groups[0][n1-1]+1
                            i = pairmap[j]
                    else:
                        n1 = len(loop_groups[0])
                        n2 = len(loop_groups[1])

                        i = loop_groups[0][0]-1
                        j = pairmap[i]

                        p = loop_groups[0][n1-1]+1
                        q = pairmap[p]

                    type1 = vienna_parameters.pair_type(sequence[i],sequence[j])
                    type2 = vienna_parameters.pair_type(sequence[q],sequence[p])

                    elements[ii].score_ = vienna_parameters.loop_energy(n1,n2,type1,type2,
                            vienna_parameters.letter_to_sequence_type(sequence[i+1]),
                            vienna_parameters.letter_to_sequence_type(sequence[j-1]),
                            vienna_parameters.letter_to_sequence_type(sequence[p-1]),
                            vienna_parameters.letter_to_sequence_type(sequence[q+1]), True, True) / 100.0
                else:
                    if(num_loop_groups > 0):
                        elements[ii].score_ = vienna_parameters.ml_energy(pairmap,sequence,loop_groups[0][0],False) / 100.0
                    else:
                        parent_elem = elements[ii].parent_

                        if(parent_elem.type_ != RNAELEMENT_STACK):
                            print "ERROR: Multiloop parent is not a stack"
                            sys.exit(0)

                        elements[ii].score_ = vienna_parameters.ml_energy(pairmap,sequence,parent_elem.indices_[len(parent_elem.indices_) -2] + 1,False) / 100.0

        elif(elements[ii].type_ == RNAELEMENT_STACK):
            indices = elements[ii].indices_
            elements[ii].quad_scores_ = []
            total_energy = 0
            type1 = vienna_parameters.pair_type(sequence[indices[0]],sequence[indices[1]])
            for jj in range(2,len(indices)):
                if(jj % 2 == 1):
                    continue
                type2 = vienna_parameters.pair_type(sequence[indices[jj+1]],sequence[indices[jj]])
                quad_energy = vienna_parameters.get_stack_score(type1,type2,True, True)
                total_energy += quad_energy
                type1 = vienna_parameters.reverse_pair_type[type2]
                elements[ii].quad_scores_.append(quad_energy)

            elements[ii].score_ = total_energy / 100.0


def get_dotplot(sequence):
    #os.system("echo " + sequence + " | ./vienna_windows_binaries/RNAfold.exe -p > rnafold_dump")
    if sequence.find("&") != -1:
        os.system("echo \"" + sequence + "\" |" + rnamake.settings.VIENNA_BIN + "RNAcofold -p > rnafold_dump")
    else:
        os.system("echo \"" + sequence + "\" |" + rnamake.settings.VIENNA_BIN + "RNAfold -p > rnafold_dump")

    # get info from output file
    try:
        file = open("dot.ps", "r")
    except IOError:
        print "Can't find dot.ps!"
        sys.exit()

    dotps = file.read()
    file.close()


    lines = re.findall('(\d+)\s+(\d+)\s+(\d*\.*\d*)\s+ubox',dotps)

    # create matrix containing index i, index j and pair probability
    dots = []

    for ii in range(0,len(lines)):
        dots.append([int(lines[ii][0]) - 1, int(lines[ii][1]) - 1, float(lines[ii][2])])

    return dots


def get_synthesized_designs_from_eterna_server(normalize_score = False, source_url = "local", noopt = False ):

    # JEEFIX
    source_url = "local"

    if source_url != "local":
        conn = httplib2.Http(".cache")
        data_url = source_url

        response, csv_data = conn.request(data_url, "GET")
    else:
        csv_f = open(settings.base_dir + "/no_validation_training/eterna_get_synthesized.csv")
        csv_data = csv_f.read()
        csv_f.close()

    csv_reader = csv.reader( csv_data.split("\n") )
    row_num = 0

    designs = []

    for row in csv_reader:
        if row_num == 0 :
            tokens = row
        else :
            col_ii = 0
            design = {}
            for item in row :
                design[tokens[col_ii]] = item
                col_ii +=1
            if("score" in design):
                if(noopt == False):
                    fill_design(design)
                designs.append(design)
        row_num += 1

    if(normalize_score):
        scores = []
        scoresum = 0
        for ii in range(0,len(designs)):
            scores.append(designs[ii]['score'])
            scoresum += designs[ii]['score']

        mean = scoresum
        if(len(designs) > 0):
            mean = mean / float(len(designs))

        stdev = 0
        for ii in range(0,len(designs)):
            stdev += (designs[ii]['score'] - mean) * (designs[ii]['score'] - mean)

        if(len(designs) > 0):
            stdev = stdev / float(len(designs))

        stdev = math.sqrt(stdev)

        for ii in range(0,len(designs)):
            designs[ii]['score'] = (designs[ii]['score'] - mean) / stdev
            designs[ii]['normalize_mean'] = mean
            designs[ii]['normalize_stdev'] = stdev

    return designs

def get_player_lab_puzzles_from_eterna_server():
    conn = httplib2.Http(".cache")
    data_url = "http://eterna.cmu.edu/get_player_lab_puzzles.php" + "?" + str(random.random());

    response, csv_data = conn.request(data_url, "GET")

    csv_reader = csv.reader( csv_data.split("\n") )
    row_num = 0

    puzzles = []

    for row in csv_reader:
        if row_num == 0 :
            tokens = row
        else :
            col_ii = 0
            puzzle = {}
            for item in row :
                puzzle[tokens[col_ii]] = item
                col_ii +=1

            if('title' in puzzle):
                puzzles.append(puzzle)
        row_num += 1

    return puzzles

def fill_design(design):
    print "Filling %s" % design['solnid']

    design['puznid'] = int(design['puznid'])
    design['gu'] = int(design['gu'])
    design['gc'] = int(design['gc'])
    design['ua'] = int(design['ua'])
    design['fe'] = float(design['fe'])

    if('score' in design and len(design['score']) > 0):
        design['score'] = float(design['score'])
    design['meltpoint'] = float(design['meltpoint'])
    design['sequence'] = design['sequence'].upper()
    if(len(design['secstruct']) + 25 == len(design['sequence'])):
        design['secstruct'] = "....." + design['secstruct'] + "...................."
    design['secstruct_elements'] = get_rna_elements_from_secstruct(design['secstruct'])
    design['pairmap'] = get_pairmap_from_secstruct(design['secstruct'])
    design['dotplot'] = get_dotplot(design['sequence'])

    fill_energy(design['secstruct_elements'],design['sequence'],design['pairmap'])
    elements = design['secstruct_elements']
    total_energy = 0.0
    for kk in range(0,len(elements)):
        total_energy += elements[kk].score_

    if(abs(total_energy - design['fe']) > 0.05):
        print "ERROR : Calculated energy doesn't match server value " + design['puztitle'] + " " + design['soltitle']
        print "SERVER : " + str(design['fe'])
        print "CALCULATED : " + str(total_energy)


def get_design_from_sequence(sequence,secstruct):
    """
    get the design dict of the sequence

    args:
    sequence is the sequence string
    secstruct is the secondary structure

    returns:
    dictionary containing various features of the design
    """
    # add basic characterstics of the design
    design = {}
    design['sequence'] = sequence
    design['secstruct'] = secstruct
    pairmap = get_pairmap_from_secstruct(secstruct)
    design['pairmap'] = pairmap

    # count pair types
    gc_count = 0
    gu_count = 0
    au_count = 0

    for ii in range(0,len(pairmap)):
        if(pairmap[ii] > ii):
            pair = sequence[ii] + sequence[pairmap[ii]]
            pair = pair.upper()

            if(pair == "GC" or pair == "CG"):
                gc_count += 1
            elif(pair == "GU" or pair == "UG"):
                gu_count += 1
            elif(pair == "AU" or pair == "UA"):
                au_count += 1
            else:
                print("Wrong pair type " + pair)

    design['gc'] = gc_count
    design['gu'] = gu_count
    design['ua'] = au_count
    design['secstruct_elements'] = get_rna_elements_from_secstruct(design['secstruct'])

    design['dotplot'] = get_dotplot(design['sequence'])
    fill_energy(design['secstruct_elements'],design['sequence'],design['pairmap'])
    elements = design['secstruct_elements']
    total_energy = 0.0
    for kk in range(0,len(elements)):
        total_energy += elements[kk].score_
    design['fe'] = total_energy
    design['meltpoint'] = 97
    return design


def load_strategy_from_file(filepath):
    class_inst = None
    expected_class = 'Strategy'

    mod_name,file_ext = os.path.splitext(os.path.split(filepath)[-1])

    if file_ext.lower() == '.py':
        py_mod = imp.load_source(mod_name, filepath)

    elif file_ext.lower() == '.pyc':
        py_mod = imp.load_compiled(mod_name, filepath)

    if expected_class in dir(py_mod):
        class_inst = py_mod.Strategy()

    return class_inst

