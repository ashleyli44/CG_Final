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

# Quickstart Instructions
1. Install mmh3, bitarray, and biopython
```
pip install mmh3
pip install bitarray
pip install biopython
```

2. Add reference fasta files to REFERENCE_FOLDER, and add contaminant fasta files to CONTAMINANT_FOLDER

3. Run main
```
python main.py data/contaminants/ data/ref_genome/ data/test/
```
To adjust any parameters of the sequence bloom tree change the desired parameter in the call to the main method at the bottom of main.py. To adjust what dataset is queried, change the final parameter in the call to the main method (True for quering sequencing data and false for querying genomic data). 

# Installation Issues
Problem: Installing dependencies fails on Windows machines
- Solution: Download Visual C++ Build Tools from https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2019
