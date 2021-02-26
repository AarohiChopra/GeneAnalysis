# gene_analysis.py
# Author: Quynh Le and Aarohi Chopra
# Updated:December3,2018
# Purpose:To analyse the hemoglobin beta gene
# Program uses function...

import re

# Part 4)
file = open("HBB_gene.txt",'r+') # opening the file HBB_gene.txt 
HBB_sequence = file.read() # reading the file
file.close() # closing file

def get_base_percentage(dna): # creating a function get_base_percentage
    dna = dna.upper() # converting all the characters in upper string
    length = len(dna) # calculating length of the variable dna
    a_content = dna.count("A") # counting the number of A's
    a_percent = round((a_content/length)*100,2) # calculating the percentage of "A"s in the dna and rounding it to 2 decimal places
    t_content = dna.count("T") # counting the number of T's
    t_percent = round((t_content/length)*100,2) # calculating the percentage of "T"s in the dna and rounding it to 2 decimal places
    c_content = dna.count("C") # counting the number of C's
    c_percent = round((c_content/length)*100,2) # calculating the percentage of "C"s in the dna and rounding it to 2 decimal places
    g_content = dna.count("G") # counting the number of G's
    g_percent = round((g_content/length)*100,2) # calculating the percentage of "G"s in the dna and rounding it to 2 decimal places
    return(str(a_percent) + "% of A's, " + str(t_percent) + "% of T's, " + str(c_percent) + "% of C's and " + str(g_percent) + "% of G's.") # returning all the percentages calculated
print("Hemoglobin beta gene has " + str(get_base_percentage(HBB_sequence)),"\n") # print the percentages


# Part 5) Count dinucleotides
# list of all the dinucleotides
dinucleotides = ['AA','AT','AG','AC',
                 'TA','TT','TG','TC',
                 'GA','GT','GG','GC',
                 'CA','CT','CG','CT']
# creating an emtpy list
all_counts = {}
for dinucleotide in dinucleotides: # starting a for loop
    count = HBB_sequence.count(dinucleotide) # counting the number of dinucleotides in the HBB_sequence 
    # We create a dictionary where dinucleotide is the key
    # and count is the value
    all_counts[dinucleotide] = count
# The dictionary is a big improvement over the lists (previously seen)
# dinucleotides and their counts are stored together as pairs key:value
print(all_counts,"\n") # printing all counts

# Part 6) Count trinucleotides
bases = ['A','T','G','C'] # creating  a list of all the trinucleotides
tri_counts = {} # creating an empty dictionary
for base1 in bases: # creating nested loops 
    for base2 in bases:
        for base3 in bases:
        # create a Trinucleotide
            trinucleotide = base1 + base2 + base3 # adding the elements to create the trinucleotides 
            count = HBB_sequence.count(trinucleotide) # counting the number of trinuleotides in the file "HBB_gene.txt"
            if count > 0: # using conditional statement for a specific condition
                tri_counts[trinucleotide] = count #adding it to the dictionary 
print(tri_counts,'\n') # printing the tricounts


# Part 7) Get AT content of HBB gene
def get_at_content(dna): # creating a function to calculate the AT contents
    """This function returns the content of AT in a DNA sequence""" 
    dna = dna.upper() # converting the DNA into upper case
    a_content = dna.count("A") # counting the number of A's
    t_content = dna.count("T") # counting the number of T's
    at_total = a_content + t_content # calculating the total number of A's and T's 
    length = len(dna) # calculating the length of the DNA
    at_content = (at_total/length) # AT contents
    return round(at_content,2) # returning the AT content upto 2 decimal places
print("The AT content of hemoglobin beta gene is " + str(get_at_content(HBB_sequence)) + ".\n") # calling the function to print AT contents of the file "HBB_sequence.txt"


# Part 8)
# Write exons into an output file called coding_sequence.txt 
exons_file = open("exons.txt",'r') # opening a file called "exons.txt" 
exons_positions = exons_file.read() # reading the file
exons_list = exons_positions.split('\n') # creating a list of all the lines as elements 
exons_file.close() # closing the file

HBB_coding_sequence = ""  # declaring an empty string as HBB_coding_sequence 
for exons in exons_list: # staring a for loop
    if(exons !=""):
        start = int(exons.split(',')[0]) #first spliting into a list then using the first element of the list as the starting value
        stop = int(exons.split(',')[1])  #spliting into a list and then using the second element as the ending value
        exon = HBB_sequence[start-1:stop] #exracting elements using the "start" and "stop" as the indexes
        HBB_coding_sequence = HBB_coding_sequence + exon # adding exons to HBB_coding_sequence

output = open("HBB_coding_sequence.txt", "w") # opening a file HBB_coding_sequence.txt
output.write(HBB_coding_sequence) # opening the file in write mode
output.close() # closing the file
print("The coding sequence of hemoglobin beta gene is: " + HBB_coding_sequence + ".\n") # printing the HBB coding sequence for dna

# Write introns into an output called nonconding_sequence.txt
introns_file = open("introns.txt",'r') # opening a file introns.txt 
introns_positions = introns_file.read() # reading the file
introns_list = introns_positions.split('\n') # creating a list called introns_list which contains all the lines in introns_position as every element
introns_file.close() # closing the file

HBB_noncoding_sequence = "" # declaring an empty string as HBB_noncoding_sequence
for introns in introns_list: # starting a for loop to iterate over the list
    start = int(introns.split(',')[0]) # start is the first element in the newly created list
    stop = int(introns.split(',')[1]) # stop is the second element in the created list
    intron = HBB_sequence[start-1:stop] # the elements between the start and stop position
    HBB_noncoding_sequence = HBB_noncoding_sequence + intron # adding the elements 

output = open("HBB_noncoding_sequence.txt", "w") # opening a file "HBB_noncoding_sequence.txt" 
output.write(HBB_noncoding_sequence.lower()) # writing the non coding sequence in the file HBB non coding sequence 
output.close() # closing the file
print("The noncoding sequence of hemoglobin beta gene is: " + HBB_noncoding_sequence.lower() + ".\n") # printing the non coding sequence


# Find all starting positions (ATG)
starting_codons = re.finditer(r"ATG",HBB_coding_sequence) # creating a list using reqular expression to find "ATG" in the coding sequence
ATG_list = [] # creating an empty list ATG_list
for codons in starting_codons: # staritng a loop 
    starting_positions = codons.start()  # finding the starting position of codons
    ATG_list.append(starting_positions) # adding the starting positions in the ATG list
print("The starting positions of the coding sequence of hemoglobin beta gene are:",ATG_list,"\n") # printing the starting position

# Translate HBB coding sequences
def translate(seq): # creating a variable
    """Translate a DNA sequence to an amino acid sequence."""
    # creating a dictionary of all the geneticode
    geneticode = { 
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    length = len(seq) # calulating the length of the sequence
    pos = 0 # initialising the variable pos as 0
    # Let us save the amino acid sequence in a list called protein
    protein = [] # creating an empty list 
    while (pos < (length-2)):       # while we still have codons to translate
        codon = seq[pos:pos+3] # Taking the first three letters 
        # Get the appropriate amino acid from the dictionary
        aa = geneticode[codon] # matching the value with
        protein.append(aa) # adding aa into the list "protien"
        pos = pos + 3 # increasing the value of the position by 3
    return "".join(protein) # returning protien
print("Translation of the hemoglobin gene coding sequence is: " + translate(HBB_coding_sequence) + ".\n")


# Where is the mutation in sickle cell anemia?
sickle_cell_file = open("sickle_cell_anemia.txt","r") # opening a file 
sickle_cell_sequence = sickle_cell_file.read() # reading the file
def unique_character(gene_A, gene_B): # creating a function
    length_A = len(gene_A) # calculating the length of gene A
    length_B = len(gene_B) # calculating the length of gene B
    mutation_position = [] # creating an empty list
    end = min(length_A,length_B) # end is the variable that stores the lesser length 
    for p in range(0, end): # creating the loop in which iterates 
        if gene_A[p] != gene_B[p]: # if the 
            mutation_position.append(p) # adding the unmatched position to the empty list
    return mutation_position # returning mutation_position
print("The position(s) of nucleotide(s) in which the mutated sickle cell gene and normal hemoglobin beta gene differ is/are:",\
      unique_character(HBB_coding_sequence,sickle_cell_sequence)) # printing The positions of nucleotides in which the mutated sickle cell gene and normal hemoglobin beta gene differ
print("The translation of sickle cell gene is: " + translate(sickle_cell_sequence) + ".\n") # printing the translation of sickle cell gene

# Where is the mutation in beta thalasemia?
beta_thalasemia_file = open("beta_thalasemia.txt","r") # opening a file
beta_thalasemia_sequence = beta_thalasemia_file.read() # reading the file
print("The position(s) of nucleotide(s) in which the mutated beta-thalasemia gene and normal hemoglobin beta gene differ is/are:", \
      unique_character(HBB_coding_sequence,beta_thalasemia_sequence)) # printing The positions of nucleotides in which the mutated beta-thalasemia gene and normal hemoglobin beta gene differ 
print("The translation of beta thalasemia gene is: " + translate(beta_thalasemia_sequence) + ".\n") # printing the translation of beta-thalasemia gene

# What happen if G in position 5 of the coding sequence is changed to A?

HBB_list = list(HBB_coding_sequence) # converting in list
HBB_list[5] = "A" # Assigning the fifth index of HBB_list to "A"
HBB_mutated_sequence = "".join(HBB_list) # Converting HBB_list into a string using the join function
print("The position(s) of nucleotide(s) in which the mutated HBB gene and normal hemoglobin beta gene differ is/are:",\
      unique_character(HBB_coding_sequence,HBB_mutated_sequence)) # printing The positions of nucleotides in which the mutated HBB gene and normal hemoglobin beta gene differ
print("The position(s) of amino acid in which the mutated HBB protein above and normal hemoglobin beta protein differ is/are:",\
      unique_character(translate(HBB_coding_sequence),translate(HBB_mutated_sequence)),"\n") # printing The positions of amino acid in which the mutated HBB protein above and normal hemoglobin beta protein differ



