"""
Bloom filter class
"""

import mmh3
from bitarray import bitarray

class BloomFilter:
    def __init__(self, str_id, kmer_list, bf_size, h):
        self.name = str_id  # string name
        self.kmer_count = len(kmer_list)
        self.BFsize = bf_size 
        self.BF = bitarray(self.BFsize) # bitarray vector that represents the bloom filter
        self.BF.setall(False)
          
        self.h = h  # number of hash functions to use per entry
        
        for kmer in kmer_list:
            self.add(kmer)
    
    def compareTo(self, bf):
        """
            Calculates the hamming distance b/w itself and another bf
            @param bf is the bloom filter object to pass in --> assume the length is the same
            @return is the integer representing the hamming distance
            w/ 0 being the minimum and len(bf) being the maximum
        """
        # exclusive OR and then calculates the number of set bits
        return (self.BF ^ bf.BF).count(True)

    def add(self, kmer):
        for i in range(self.h):
            hashValue = mmh3.hash(kmer, i) % self.BF.length()
            self.BF[hashValue] = True

    def contains(self, kmer):
        for i in range(self.h):
            hashValue = mmh3.hash(kmer, i) % self.BF.length()
            if not self.BF[hashValue]:
                return False
        return True

    def union(self, other_bloom_filter):
        if self.name == "":
            self.name = other_bloom_filter.name
        else:
            self.name = self.name + "U" + other_bloom_filter.name
        for i in range(len(self.BF)):
            self.BF[i] = self.BF[i] or other_bloom_filter.BF[i]
            
    def getBF(self):
        return self.BF
    
    def getName(self):
        return self.name    
