# CG_Final
Sequence Bloom Trees for Contaminant Detection

Final Project for 601.447 Computational Genomics: Sequences (Fall 2019)

#### Group 20 Team Members
Eric Cao, Sandeep Kambhampati, Ashley Li, Robert Li 

# Dependencies
- Python3
- MurmurHash3
- Bitarray
- Biopython

# Installation Instructions
1. Install mmh3, bitarray, and biopython
```
pip install mmh3
pip install bitarray
pip install biopython
```

2. Put fasta files in FOLDER_NAME
   Need separate folders for the reference genome and contaminant genome

3. Run main (command-line command: python main.py data/contaminants/ data/ref_genome/ data/test/)
