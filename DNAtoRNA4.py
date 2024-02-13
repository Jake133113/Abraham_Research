#Python script to get convert the DNA sequences in RefSeq to RNA Sequences. In this case, we are also give the 
#amino acid sequence to verify the length 

import re

#We want to input the sequence of letters, and get rid of all the spaces and numbers

test_sequence = input("Copy and Paste the full DNA sequence from RefSeq Database: ")
protein_sequence = input("Copy and Paste the full protein sequence from RefSeq Database: ")
protein_sequence = str(protein_sequence)

def DNAtoRNA(DNA):
        DNA = DNA.replace(' ','')
        #Now, we want to remove all numbers from the sequence
        pattern = r'[0-9]'
        DNA = re.sub(pattern, '', DNA)
        DNA = DNA.upper()
        DNA = DNA.replace('T','U')
        return str(DNA)


#Now we want to pick out the index of the first start codon (AUG), and delete it and everything before$

RNA_seq = str(DNAtoRNA(test_sequence))
index = RNA_seq.index('AUG')
index = int(index)

RNA_seq_nostart = RNA_seq[index+3:]
#print(RNA_seq_nostart)

#creating a for loop that scans every third index for one of the stop codons.
def Stop_sequence(sequence):
    j = 0
    while j<len(sequence):
        if sequence[j:j+3] != 'UAA' and sequence[j:j+3] != 'UAG' and sequence[j:j+3] !='UGA':
            j+=3
        elif sequence[j:j+3] == 'UAA':
            final_sequence = sequence[:j]
            break
        elif sequence[j:j+3] == 'UAG':
            final_sequence = sequence[:j]
            break
        elif sequence[j:j+3] == 'UGA':
            final_sequence = sequence[:j]
            break
    return final_sequence


test1 = Stop_sequence(RNA_seq_nostart)
    
if len(test1) == 3*len(protein_sequence):
    print("The coding sequence is the correct length for the resulting protein sequence")
    print("Final sequence is: ", test1)
else:
    print('The sequence was the incorrect length. This code will now try adding the start codon right before the first detection of the start codon because the start codon also codes for the amino acid "M". If that does not work, it will try adding the start codon to the very beginning.')

start_codon = 'AUG'
test2 = Stop_sequence('AUG' + RNA_seq_nostart)
test3 = Stop_sequence(RNA_seq)

if len(test2) == 3*len(protein_sequence):
        print('The sequence is now the correct length with a start codon added right before the first occurance of "AUG". ')
        print('Final sequence is: ', test2)
elif len(test3) == 3*len(protein_sequence):
        print('The sequence is now the correct length with a start codon added right in the beginning. ')
        print('Final sequence is: ', test3)
else:
        print('The codon sequence is still the wrong length.') 


