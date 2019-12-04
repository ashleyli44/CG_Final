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
        @return a test dictionary {species_name_num : read}
    '''
    species = list(data_dict) #keys of the dictionary = species
    ret_dict = {}
    
    num_entries = int(n/len(species)) #number of entries per species
    #array to hold the sequences
    #each row is a diff species
    #each col is a diff seq for that species
    species_seq = np.zeros((len(species), num_entries ), dtype = 'U1000000')

    #iterate through the data_dict by species
    for i in range(0, len(species)):
        reads = [] 
        sub_dict = data_dict[species[i]]
        if (len(list(sub_dict))) < num_entries: #tests if we need to split the reads
            str_all = sub_dict.values() #all the values
            n_split = int(num_entries / len(list(sub_dict))) #how many times I need to split each read
            if n_split == 1: #avoids situation where n_split is rounded too low
                n_split += 1
            for j in str_all: #iterate through the reads to split them
                len_split = int(len(j)/n_split) #how long each split should be based on the number of splits
                reads.append([(j[x:x+len_split]) for x in range(0, len(j), len_split)])
            reads = [item for sublist in reads for item in sublist] #flattens list
        else: #if there are enough reads to randomly sample, then do so
            idx_keys = random.sample(list(sub_dict), int(n/len(species))) 
            reads = [sub_dict[x] for x in idx_keys]
        
        species_seq[i,:] = np.array(reads[0:num_entries], dtype = str) 
        #reads may be longer than n_entries depending on how you split it

    #creation of the return dictionary as dictionary of dictionaries, same as good data set            
    ret_dict = {} 
    for i in range(0, len(species)):
      for j in range(0, num_entries):
            idx = species[i] + "_" + str(j)
            ret_dict[idx] = species_seq[i,j]
    
    return(ret_dict)


bac_folder = sys.argv[1]
ref_folder = sys.argv[2]

#reads in fasta files
bac_dict = dataParse(bac_folder)
ref_dict = dataParse(ref_folder)

#creates datasets
n = 500
good_dict = goodData(ref_dict, n)
bad_dict = badData(bac_dict, n)

