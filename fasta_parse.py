from Bio import SeqIO
import sys
import os
import numpy as np
import random

def fastaParse(input_file):
    '''
        Function parses a fasta file
        @param input_file is the filename
        @return a dictionary where the keys are seq IDs and the values are sequences
    '''
    fasta_sequences = SeqIO.parse(open(input_file),'fasta')
    fasta_dict = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        fasta_dict[name] = sequence
    return(fasta_dict)

def dataParse(input_dir):
    '''
        Function automatically parses all the files in a data folder
        @param input_dir is the name of the directory to extract all the files
        @return a dictionary of dictionarys {species_name: {seq_ID: seq}, ... }
    '''
    data_dict = {}
    for filename in os.listdir(input_dir):
        ref_id = os.path.splitext(filename)[0]
        ifile = input_dir + filename
        data_dict[ref_id] = fastaParse(ifile)
    return(data_dict)

def goodData(data_dict, n):
    '''
        Function creates the good data set from the reference genome
        @param data dict is the input dictionary
        @param n is the number of reads to generate (must be divisible by 5)
        @return a test dictionary {species_name_num : read}
        Note: reads can be of variable length
    '''
    species = list(data_dict) #species in the reference genome (just one in this case)
    keys = random.sample(list(data_dict[species[0]]), n) #randomly samples keys
    ret_dict = {}
    for i in range(0, n):
        key = species[0] + "_" + str(i) #creates a key name ex. mouse_1 --> <species>_<i>
        ret_dict[key] = data_dict[species[0]][keys[i]]
    return(ret_dict)

def badData(data_dict, n):
    '''
        Function creates the contaminant set
        @param data_dict is the input dictionary
        @param n is the number of reads to generate (must be divisible by 5)

    '''
    species = list(data_dict)
    m = len(species)/n
    ret_dict = {}
    
    num_entries = int(n/len(species))
    keys = np.zeros((len(species), num_entries ), dtype = 'U1000000')
    for i in range(0, len(species)):
        reads = []
        sub_dict = data_dict[species[i]]
        if (len(list(sub_dict))) < num_entries:
            str_all = sub_dict.values()
            n_split = int(num_entries / len(list(sub_dict))) #how many times I need to split each read
            if n_split == 1:
                n_split += 1
            for j in str_all:
                len_split = int(len(j)/n_split) 
                reads.append([(j[x:x+len_split]) for x in range(0, len(j), len_split)])
            reads = [item for sublist in reads for item in sublist]
        else:
            idx_keys = random.sample(list(sub_dict), int(n/len(species)))
            reads = [sub_dict[x] for x in idx_keys]
        keys[i,:] = np.array(reads[0:num_entries], dtype = str)
        keys[i,0] = reads[0]

    for i in range(0, len(species)):
        for j in range(0, num_entries):
            idx = species[i] + "_" + str(j)
            ret_dict[idx] = keys[i,j] 
    return(ret_dict)

bac_folder = sys.argv[1]
ref_folder = sys.argv[2]

bac_dict = dataParse(bac_folder)
ref_dict = dataParse(ref_folder)

n = 50

good_dict = goodData(ref_dict, n)
bad_dict = badData(bac_dict, n)

