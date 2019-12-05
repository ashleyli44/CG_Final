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
import math
import numpy as np
import random
from Bio import SeqIO

def main(myKmer_size, myFpr, myThreshold):

    bac_folder = sys.argv[1]
    ref_folder = sys.argv[2]

    #reads in fasta files
    bac_dict = dataParse(bac_folder)
    ref_dict = dataParse(ref_folder)
    
    #creates datasets
    n = 500
    good_dict = goodData(ref_dict, n)
    bad_dict = badData(bac_dict, n)

    kmer_size = myKmer_size	
    def convertReadtoKmerList(read):
      kmer_list = []
      for i in range(len(read)-kmer_size + 1):
        kmer = read[i:i+kmer_size]
        kmer_list.append(kmer)
      return kmer_list

    total_kmerList = []
    for species in bac_dict:
      species_dict = bac_dict[species]
      for readKey in species_dict:
        read = species_dict[readKey]
        kmerList = convertReadtoKmerList(read)
        total_kmerList = total_kmerList + kmerList

    fpr = myFpr
    BFsize = getBFsize(len(total_kmerList), fpr)
    BFHashCount = getHashFunctionCount(len(total_kmerList), BFsize)
    #print("BFsize: " + str(BFsize))
    #print("Hash Count: " + str(BFHashCount))
    
    print("Constructing species BFs")
    BFList = []
    for species in bac_dict:
      print(species) 
      species_dict = bac_dict[species]
      kmerList = []
      for readKey in species_dict:
        #print(readKey)
        read = species_dict[readKey]
        kmers = convertReadtoKmerList(read)
        kmerList = kmerList + kmers
      species = str(species)
      addBF = BloomFilter(species, kmerList, BFsize, BFHashCount) 
      BFList.append(addBF) 

    print("Constructing Bloom Tree")

    from bloom_tree import BloomTree
    from bloom_node import BloomNode
    inverseBloomTree = BloomTree(myThreshold, BFsize, BFHashCount)
    for bloomFilter in BFList:
      newNode = BloomNode(bloomFilter)
      inverseBloomTree.add(newNode)

    bloomTreeSize = sys.getsizeof(inverseBloomTree)
    numTotal = 0
    numAccurate = 0
    totalMatches = 0
    print("Querying")
    #print(bad_dict.keys())
    for species in bad_dict:
      species_BadDict = bad_dict[species]
      for readID in species_BadDict:
        read = species_BadDict[readID]
        queryKmers = convertReadtoKmerList(read)
        possibleMatches = inverseBloomTree.query(queryKmers)          
        possibleMatchNames = []
        for match in possibleMatches:
          #print(match.bloom_filter)
          possibleMatchNames.append(match.bloom_filter.getName())
        speciesName = str(species)
        if speciesName in possibleMatchNames:
          numAccurate += 1
          numTotal += 1
          totalMatches += len(possibleMatchNames)
        else:
          numTotal += 1
          totalMatches += len(possibleMatchNames)
    accuracy = numAccurate/numTotal
    avgMatches = totalMatches/numTotal
    print("Accuracy: " + str(accuracy))
    print("AvgMatches: " + str(avgMatches))        
   
    numMouseMatches = 0
    numQueries = 0
    for readKey in good_dict:
      read = good_dict[readKey]
      queryKmers = convertReadtoKmerList(read)
      possibleMatches = inverseBloomTree.query(queryKmers)
      numMouseMatches += len(possibleMatches)
      numQueries += 1
    avgMouseMatches = numMouseMatches/numQueries
    print("avgMouseMatches: " + str(avgMouseMatches))  
   
    print("BF size: " + str(BFsize))
    print("Hash Count: " + str(BFHashCount))
    print("Bloom Tree size: " + str(bloomTreeSize))     

def getBFsize(n,fpr):
    m = -(n * math.log(fpr)) / (math.log(2) ** 2)
    return int(m)

def getHashFunctionCount(n, m):
    k = (m / n) * math.log(2)
    return int(k)

main(5, 0.1, 0.3)
