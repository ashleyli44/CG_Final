from Bio import SeqIO
import sys

input_file = sys.argv[1]

def fastaParse(input_file):
    fasta_sequences = SeqIO.parse(open(input_file),'fasta')
    fasta_dict = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        fasta_dict[name] = sequence
    return(fasta_dict)

print(fastaParse(input_file))
