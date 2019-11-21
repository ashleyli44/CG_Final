from Bio import SeqIO
import sys
import os


def fastaParse(input_file):
    fasta_sequences = SeqIO.parse(open(input_file),'fasta')
    fasta_dict = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        fasta_dict[name] = sequence
    return(fasta_dict)

data_folder = sys.argv[1]

data_dict = {}

for filename in os.listdir(data_folder):
    ref_id = os.path.splitext(filename)[0]
    ifile = data_folder + filename
    data_dict[ref_id] = fastaParse(ifile)





