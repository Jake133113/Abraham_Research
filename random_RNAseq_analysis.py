#Python script that generates random RNA sequences at different lengths from 50 to 3000 bases and extracts all possible stem formations in each. 
#At each length, the average stem length is found along with its SEM and the average maximum stem length is found with its SEM.

import pandas as pd
from scipy import stats
import numpy as np

bases = ['C','G','A','U']

#Initialize empty arrays that will be filled with random RNA sequences at each specified length
L50 = []
L100 = []
L150 = []
L200 = []
L250 = []
L300 = []
L500 = []
L1000 = []
L2000 = []
L3000 = []

def gen_rand_sequence(len):
        #this function generates a random RNA sequence at a specified length (in number of bases)
        seq = ''
        for i in range(len):
                num = np.random.randint(0,4)
                base = bases[num]
                seq = seq + base
        return seq

def gen_seq_list(num_seqs, seq_len, list_name):
        #this function generates a list of random sequences. Inputs give the number of sequences, how long they are, and the list name. 
        for i in range(num_seqs):
                rand_seq = gen_rand_sequence(seq_len)
                list_name.append(rand_seq)
        return list_name

#specify desired sequence lengths
seq_lengths = [50,100,150,200,250,300,500,1000,2000,3000]
sequence_lists = [L50,L100,L150,L200,L250,L300,L500,L1000,L2000,L3000]

for i in range(len(seq_lengths)):
        gen_seq_list(100,seq_lengths[i],sequence_lists[i])

def _calc_stem_separation(first_base, last_base, stem_length):
        #Function adapted from Dillion Fox and Brian Andrews Quvax Repository for RNA folding optimization
        return (last_base - (stem_length - 1)) - (first_base + (stem_length - 1)) - 1

WC_interactions = [("A", "U"),("U", "A"),("G", "C"),("C", "G")]
WC_wobble_interactions = [("A", "U"),("U", "A"),("G", "C"),("C", "G"),("G", "U"),("U", "G")]

def gen_stems(seq_len, min_loop, min_stem, span, sequence, interactions):
        #Function adapted from Dillion Fox and Brian Andrews Quvax Repository for RNA folding optimization
        #Inputs are sequence length, minimum loop length (recommended is 3), minimum stem length (recommended is 3), 
        #span (sets the distance of which bases can form bonds, 0 for no span), the RNA sequence,
        #and interactions (can include wobble or only Watson-Crick interactions)
        """
        Generates a list of tuples with three entries:
        - first two elements are sequence indices of a base pair identifying the stem
        - third element is the length of the stem
        - Ex: stems = [(1, 13, 3), ()]
            stems[0] = (1, 13, 3) corresponds to stem with base pairs between
            nucleotides with index 1 and 13, 2 and 12, 3 and 11.

        """
        pairs = []
        for i in range(seq_len - 2 * min_stem - min_loop):
            for j in range(i + 2 * min_stem + min_loop - 1,seq_len):
                for k in range(1, seq_len):  # start at 1 because can't have stem of length 0
                    if i + k >= seq_len:
                        break
                    if _calc_stem_separation(i, j, k) < min_loop or (span > 0
                        and _calc_stem_separation(i, j, k) > span):
                        break
                    # (k-1) in following conditional because stem (1,5,2) has pairs (1,5),
                    # (2,4). i+k = 1+2 = 3 and base 3 is not in the stem
                    if (sequence[i + (k - 1)],sequence[j - (k - 1)]) in interactions:
                        if k >= min_stem:
                            # increment i and j because they start index at zero
                            # and RNA structure files index at 1
                            # don't increment stem length because leads to
                            # overriding min loop length requirement in special cases
                            pairs.append((i + 1, j + 1, k))
                    else:
                        break
        stems = pairs
        return stems

def generate_stem_data(seq_data):
    max_stems = []
    avg_stems = []
    print("Analyzing New Set of Sequences")
    count = 0
    for sequence in seq_data:
        count += 1
        print("Analyzing Sequence: ", count)
        lengths = []
        seq = sequence
        stems = gen_stems(len(seq),3,3,0,seq,WC_wobble_interactions)
        max_stems.append(max(stem[2] for stem in stems))
        for stem in stems:
            lengths.append(stem[2])
        avg_stems.append(stats.tmean(lengths))
    #return avg stem length, standard error of that mean, avg max stem length, standard error of the max 
    avg_stem = stats.tmean(avg_stems)
    avg_stem_error = stats.sem(avg_stems)
    avg_max_stem = stats.tmean(max_stems)
    max_stem_error = stats.sem(max_stems)
    return (avg_stem, avg_stem_error, avg_max_stem, max_stem_error)

final_data = []

for seq_data in sequence_lists:
    final_data.append(generate_stem_data(seq_data))

print("Output for these sequence lengths: ", seq_lengths)
print("For each length, the output will be: (average stem, avg stem standard error, average maximum stem, max stem standard error)")
print(final_data)

