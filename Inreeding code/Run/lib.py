#!/usr/bin/env python
from FisSimFunctions import *

def get_output_pairs(input_list):
    """
    Given an input list of format [1, 2, 3] returns a set of all the possible
    two digit combinations. The order of the combinations doesn't matter, so
    the output is a set of pairs.
    """
    output = set()
    for n1 in input_list:
        # TODO: be cleverer about this loop to avoid unnecessary computation.
        for n2 in input_list:
            # We sort the pair to make sure that the ordering doesn't matter.
            # The set won't duplicate the element if we've already inserted it
            # this way.
            pair = tuple(sorted((n1, n2)))
            output.add(pair)
    return output


def get_combinations(input_genes, num_output_pairs):
    """
    Given a list of integers `input_genes`, return every possible combination
    of the genes that fit into the given number of `num_output_pairs`. The
    output pairs are sets because the order of genes within the pairs does not
    matter. The order of the pairs *does* matter.

    Notes:
    - Each gene in `input_genes` must be present at least one output pair.
    """
    # Input validation
    if (num_output_pairs * 2) < len(input_genes):
        raise ValueError("num_output_pairs x 2 must be greater than the number of input genes")

    possible_output_pairs = list(get_output_pairs(input_genes))

    # We need to find all possible combinations of the output pairs
    # that fit into num_output_pairs
    output = []

    # General algorithm
    # Step 1: make literally every combo
    # Step 2: prune ones that don't have every digit represented

    # Think of this list like a combination lock. We're going to increment each
    # index into `possible_output_pairs` one by one until we've gone through
    # every possible combination.
    pair_indexes = [0] * num_output_pairs

    # while all indexes not == (to num_output_pairs - 1)
    num_pairs = len(possible_output_pairs)
    # when we're done incrementing all of the indexes the final time, the first
    # index will be out of bounds of the array so we stop then
    while pair_indexes[0] < num_pairs:
        # print(pair_indexes)  # DEBUG
        # Create the list with the current values and append to output
        pair_list = []
        for index in pair_indexes:
            pair_list += [possible_output_pairs[index]]
        output += [pair_list]

        # Increment the next pair
        increment_index(pair_indexes, num_pairs)

    # PRUNE BABY
    real_output = []
    for pair_list in output:
        if contains_all_genes(pair_list, input_genes):
            real_output += [pair_list]

    return real_output


def contains_all_genes(pair_list, input_genes):
    # THIS IS HORRIBLY INEFFICIENT - fix me first if you need speed
    for num in input_genes:
        found = False
        for pair in pair_list:
            if num in pair:
                found = True
                break
        if not found:
            return False
    return True



def increment_index(pair_indexes, num_pairs):
    """
    Given a list of pair indexes and the number of pairs to choose from,
    increments the right most index by one, and all indexes to the left of
    that as appropriate from the overflow.

    Mutates the list given and also returns it.
    """
    num_output_pairs = len(pair_indexes)

    index = 1
    while True:
        pair_indexes[num_output_pairs - index] = pair_indexes[num_output_pairs - index] + 1
        # Don't go to index -1
        if index >= num_output_pairs:
            break
        # If the index is still a valid pair-list index and we don't overflow, return
        if pair_indexes[num_output_pairs - index] < num_pairs:
            break
        # We've gone over the max of the pairs list, so set to 0 and go one index to the left to increment that next
        pair_indexes[num_output_pairs - index] = 0
        index += 1

    return pair_indexes


def increment_index2(pair_indexes, max_per_index):
    """
    Same as above but allows different maxes for each index.
    """
    num_output_pairs = len(pair_indexes)

    index = 1
    while True:
        pair_indexes[num_output_pairs - index] = pair_indexes[num_output_pairs - index] + 1
        # Don't go to index -1
        if index >= num_output_pairs:
            break
        # If the index is still a valid pair-list index and we don't overflow, return
        if pair_indexes[num_output_pairs - index] < max_per_index[num_output_pairs - index]:
            break
        # We've gone over the max of the pairs list, so set to 0 and go one index to the left to increment that next
        pair_indexes[num_output_pairs - index] = 0
        index += 1

    return pair_indexes


def tuple_same(in_tuple):
    if in_tuple[0] == in_tuple[1]:
        return 0
    return 1

# Returns a list where the value is 0 for each index where the pair has different values
# and the value is 1 if the pair is the same, for each pair in the list
# Takes a list of lists and operates on those lists
def tuples_same(pair_list_list):
    return [[tuple_same(pair) for pair in pair_list] for pair_list in pair_list_list]


# Calculate observed heterozygocity
def calculate_total_Ho(subject_h_list):
    # Input is a list of lists of lists
    # Algorithm:
    # Take first index of first subject and compare to first index of all other subjects, count number of 1's vs total
    # Take first index of first subject and compare to first index of all other subjects except add one to the final subject
    # Add two to the final
    # When wrapping, add one to the next
    # Reframe as algorithm - add 1 to 0 then spill over to the right
    num_subjects = len(subject_h_list)
    num_possibilities = []
    for l in subject_h_list:
        num_possibilities += [len(l)]
    
    subject_indexes = [0] * num_subjects
    Ho_possibilities = []
    # All subjects should be the same length so this works
    while subject_indexes[0] < num_possibilities[0]:
        # do maths
        # Collect all the numbers
        h_nums = []
        for i in range(len(subject_indexes)):
            index = subject_indexes[i]
            try:
                for h_num in subject_h_list[i][index]:
                    h_nums += [h_num]
                    #print("strooz", subject_h_list)
            except Exception as e:
                # For easier debugging
                print(num_possibilities)
                print(index)
                print(i)
                print(subject_indexes)
                print(subject_h_list[i])
                # Still raise it so it fails here
                raise e
        # Count the 1s / total
        Ho_possibilities += [h_nums.count(1) / len(h_nums)]
        
        # increment index
        increment_index2(subject_indexes, num_possibilities)

    return Ho_possibilities

def calculate_total_He(subject_h_list):
    # Input is a list of lists of lists
    # Algorithm:
    # Take first index of first subject and compare to first index of all other subjects, count number of 1's vs total
    # Take first index of first subject and compare to first index of all other subjects except add one to the final subject
    # Add two to the final
    # When wrapping, add one to the next
    # Reframe as algorithm - add 1 to 0 then spill over to the right
    num_subjects = len(subject_h_list)
    num_possibilities = []
    for l in subject_h_list:
        num_possibilities += [len(l)]
    
    subject_indexes = [0] * num_subjects
    Ho_possibilities = []
    alleleCountsList = []
    # All subjects should be the same length so this works
    while subject_indexes[0] < num_possibilities[0]:
        alleleCounts = {'total': 0}
        alleleCountsList += [alleleCounts]
        # do maths
        # Collect all the numbers
        h_nums = []
        for i in range(len(subject_indexes)):
            index = subject_indexes[i]
            try:
                for h_num in subject_h_list[i][index]:
                    h_nums += [h_num]
                    for allele in h_num:
                        if allele not in alleleCounts:
                            alleleCounts[allele] = 0
                        alleleCounts[allele] += 1
                        alleleCounts['total'] += 1
                    #print("strooz", subject_h_list)
            except Exception as e:
                # For easier debugging
                print(num_possibilities)
                print(index)
                print(i)
                print(subject_indexes)
                print(subject_h_list[i])
                # Still raise it so it fails here
                raise e
        # Count the 1s / total
        Ho_possibilities += [h_nums.count(1) / len(h_nums)]
        
        # increment index
        increment_index2(subject_indexes, num_possibilities)

    
    for alleleCounts in alleleCountsList:
        for allele, count in alleleCounts.items():
            if allele == 'total':
                continue
            alleleCounts[allele] = (alleleCounts[allele] / alleleCounts['total'])**2
            
            

    
    He_values = []
    for alleleCounts in alleleCountsList:
        value = 0
        for allele, count in alleleCounts.items():
            if allele == 'total':
                continue
            value += count
        He_values += [(1 - value)]

            
#        alleleCountsList = []
#        for possibility in range(varientNo):
#            alleleCounts = {'total': 0}
#            alleleCountsList += [alleleCounts]
#            for pair in possibility:
#                for allele in pair:
#                    if allele not in alleleCounts:
#                        alleleCounts[allele] = 0
#                    alleleCounts[allele] += 1
#                    alleleCounts['total'] += 1
    return He_values

def MakingPosibilityGraph(result, NoOfLoci, eyePrime, Fis):
    import numpy as np
    import matplotlib.pyplot as plt

    
    combinations_total = []
    tuples_same_total = []
    
    for gene_data in result:
        combinations = get_combinations(gene_data, NoOfLoci)
        tuples_same_i = tuples_same(combinations)
        print(tuples_same_i)
        combinations_total += [combinations]
        tuples_same_total += [tuples_same_i]
    
    
    
    
    ####COMBOS PER LOCI, TURING IT INTO a  LIST OF LISTS PER LOCI####
    combinations_total_oneloci = []
    combinations_total_perloci = []
    
    for i in range(NoOfLoci):
        for sublist in combinations_total:
                combinations_total_oneloci += [[[item[i]] for item in sublist]]
        combinations_total_perloci += [combinations_total_oneloci]
        combinations_total_oneloci = []            
    
    
    ####TUPS PER LOCI, TURING IT INTO a  LIST OF LISTS PER LOCI####
    tuples_total_oneloci = []
    tuples_total_perloci = []
    for i in range(NoOfLoci):
        for sublist in tuples_same_total:
                tuples_total_oneloci += [[[item[i]] for item in sublist]]
        tuples_total_perloci += [tuples_total_oneloci]
        tuples_total_oneloci = []            
    
    
    
    total_He = []
    total_Ho = []
    for i in range(NoOfLoci):
        HeHe = calculate_total_He(combinations_total_perloci[i])
        total_He += [HeHe]
        HoHo = calculate_total_Ho(tuples_total_perloci[i])
        total_Ho += [HoHo]
    
    # don't print this it'll crash
    #print(total_Ho)
    
    possible_fis_oneloci = []
    possible_fis_perloci = []
    for loci in range(NoOfLoci):
        
        for i in range(len(total_He[loci])):
            if total_He[loci][i] == 0:
                possible_fis_oneloci += [0]
            else:    
                possible_fis_oneloci += [(total_He[loci][i] - total_Ho[loci][i])/total_He[loci][i]]
        
        possible_fis_perloci += [possible_fis_oneloci]
        possible_fis_oneloci = []
    
    possible_fis = np.mean(np.array(possible_fis_perloci), axis=0)
    
        
    #with open("possible_fis.csv", "w") as f:
    #    f.write("\n".join(str(x) for x in possible_fis))
    
    #for i in range(len(tuples_same)):
    #    for l in tuples_same[i]:
    #        print(tuples_same)
    
    
    plt.hist(possible_fis,zorder=-1)
    plt.xlim((-1, 1))
    if eyePrime[0] > 1:
        eyePrime[0] = 1
    elif eyePrime[0] < -1:
        eyePrime[0] = -1
    plt.scatter(eyePrime[0],1,s = 130, c="#ffa303", marker = "^")
    plt.annotate("I'",(eyePrime[0],1), xytext = (eyePrime[0]-2,-14), textcoords = "offset pixels",zorder=1)
    plt.scatter(Fis,1, s = 130, c="#ff2503", marker = "^")
    plt.annotate("Fis",(Fis,1), xytext = (Fis-6,-14), textcoords = "offset pixels",zorder=1)
