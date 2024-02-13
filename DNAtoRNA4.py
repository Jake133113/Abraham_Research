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
            print('Final Sequence: ', sequence[:j])
            final_sequence = sequence[:j]
            break
        elif sequence[j:j+3] == 'UAG':
            print('Final Sequence: ', sequence[:j])
            final_sequence = sequence[:j]
            break
        elif sequence[j:j+3] == 'UGA':
            print('Final Sequence: ', sequence[:j])
            final_sequence = sequence[:j]
            break

print(Stop_sequence(RNA_seq_nostart))
    
#Adding in an error detection and possible correction if the length of the protein sequence isn't exactly 1/3 of the coding sequence

if len(final_sequence) == 3*len(protein_sequence):
    print("The coding sequence is the correct length for the resulting protein sequence")
else:
    print('The sequence was the incorrect length. This code will try adding the start codon to the beginning and check the length.')
    try: 
        new_seq = Stop_sequence(RNA_seq)
        print('The sequence is now the correct length with a start codon added to the beginning. ')
    except len(new_seq) != 3*len(protein_sequence):
        print('The codon sequence is still the wrong length. Please check your inputs for errors.') 
        
        
