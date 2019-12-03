"""
Bloom filter factory
"""

# READ ME: install mmh3 and bitarray via pip install mmh3 and pip install bitarray
import math
import mmh3
from bitarray import bitarray

#import numpy as np

from bloom_filter import BloomFilter


class BloomFilterFactory(object):
    
    #kmer_list should be a list, fpr should be a double
    def __init__(self, kmer_list, fpr):
        self.kmer_count = len(kmer_list)
        self.fpr = fpr
        self.BFsize = self.getBFsize(self.kmer_count, self.fpr)
        self.hashFunctionCount = self.get_hashFunctionCount(self.kmer_count, self.BFsize)
        self.BF = bitarray(self.BFsize)
        self.BF.setall(0)

        for kmer in kmer_list:
            self.add(kmer)

    @classmethod
    def getBFsize(self,n,fpr):
        m = -(n * math.log(fpr)) / (math.log(2) ** 2)
        return int(m)

    @classmethod
    def get_hashFunctionCount(self, n, m):
        k = (m / n) * math.log(2)
        return int(k)

    def add(self, kmer):
    	for i in range(self.hashFunctionCount):
          hashValue = mmh3.hash(kmer, i) % self.BFsize
          self.BF[hashValue] = True

    def checkExistence(self, item):
        for i in range(self.hashFunctionCount):
            hashValue = mmh3.hash(item, i) % self.BFsize
            if self.BF[hashValue] == False:
                return False
        return True

    def returnBF(self):
        return self.bit_array

    @staticmethod
    def new_instance(name=""):
        m = 1000
        h = 1
        return BloomFilter(name, bitarray(m), h)
