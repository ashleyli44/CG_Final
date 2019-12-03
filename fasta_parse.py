from Bio import SeqIO
import sys
import os
import numpy as np
import random
def fastaParse(input_file):
    fasta_sequences = SeqIO.parse(open(input_file),'fasta')
    fasta_dict = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        fasta_dict[name] = sequence
    return(fasta_dict)

def dataParse(input_dir):
    data_dict = {}
    for filename in os.listdir(input_dir):
        ref_id = os.path.splitext(filename)[0]
        ifile = input_dir + filename
        data_dict[ref_id] = fastaParse(ifile)
    return(data_dict)

def goodData(data_dict, n):
    species = list(data_dict)
    
    #idx = np.random.choice(list(data_dict[species]), size = n, replace = False)
    #reads = data_dict[species][idx]
    reads = random.sample(list(data_dict[species[0]]), n)
    ret_dict = {}
    for i in range(0, n):
        key = species[0] + "_" + str(n)
        ret_dict[key] = data_dict[species[0]][reads[i]]
    return(ret_dict)

bac_folder = sys.argv[1]
ref_folder = sys.argv[2]

bac_dict = dataParse(bac_folder)
ref_dict = dataParse(ref_folder)

good_dict = goodData(ref_dict, 50)
half_dict = {}
bad_dict = {}
