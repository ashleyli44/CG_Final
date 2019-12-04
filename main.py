from bloom_factory import BloomFilterFactory
from bloom_filter import BloomFilter
from bloom_node import BloomNode
from bloom_tree import BloomTree
from bitarray import bitarray

from fasta_parse import fastaParse
from fasta_parse import dataParse
from fasta_parse import goodData
from fasta_parse import badData

import sys
import os
import numpy as np
import random
from Bio import SeqIO

def main():

    bac_folder = sys.argv[1]
    ref_folder = sys.argv[2]

    #reads in fasta files
    bac_dict = dataParse(bac_folder)
    ref_dict = dataParse(ref_folder)

    #creates datasets
    n = 500
    good_dict = goodData(ref_dict, n)
    bad_dict = badData(bac_dict, n)

    kmer_size = 30	
    def convertReadtoKmerList(read):
      kmer_list = []
      for i in range(len(read)-kmer_size + 1):
        kmer = read[i:i+kmer_size]
        kmer_list.append(kmer)
      return kmer_list

    total_kmerList = []
    for species in bad_dict:
      species_dict = bad_dict[species]
      for readKey in species_dict:
        read = species_dict[readKey]
        kmerList = convertReadtoKmerList(read)
        total_kmerList = total_kmerList + kmerList

    fpr = 0.1 
    BFsize = getBFsize(length(total_kmerList), fpr)
    BFHashCount = getHashFunctionCount(length(total_kmerList),rootSize)
    
    print("Constructing species BFs")
    BFList = []
    for species in bad_dict:
      species_dict = bad_dict[species]
      kmerList = []
      for readKey in species_dict:
        read = species_dict[readKey]
        kmers = convertReadtoKmerList(read)
        kmerList = kmerList + kmers
      species = str(species)
      addBF = Bf(species, kmerList, BFsize, BFHashCount)
      BFList.append(addBF) 

    print("Constructing Bloom Tree")

    from bloom_tree import BloomTree
    from bloom_node import BloomNode
    inverseBloomTree = BloomTree(0.1)
    for bloomFilter in BFList:
      newNode = BloomNode(bloomFilter)
      inverseBloomTree.add(newNode)   

def getBFsize(n,fpr):
    m = -(n * math.log(fpr)) / (math.log(2) ** 2)
    return int(m)

def get_hashFunctionCount(n, m):
    k = (m / n) * math.log(2)
    return int(k)

main()
