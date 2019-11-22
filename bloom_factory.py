"""
Bloom filter factory
"""

# Warning: unfinalized/untested code

# READ ME: install mmh3 and bitarray via pip install mmh3 and pip install bitarray
import math
import mmh3
from bitarray import bitarray

from bloom_filter import BloomFilter


class BloomFilterFactory(object):
    def __init__(self, kmer_list, fpr):
        self.kmer_count = size(kmer_list)
        self.fpr = fpr
        self.BFsize = self.getBFsize(self.kmer_count, self.fpr)
        self.hashFunctionCount = self.get_hashFunctionCount(self.kmer_count, self.BFSize)

        self.BF = bitarray(self.m)
        self.BF.setall(0)

        for kmer in kmer_list:
            self.add(kmer)

    @classmethod
    def getBFsize(self,n,fpr):
        ##figure out math behind this
        m = -(n * math.log(fpr)) / (math.log(2) ** 2)
        return int(m)

    @classmethod
    def get_hashFunctionCount(self, m, n):
        k = (m / n) * math.log(2)
        return int(k)

    @classmethod
    def add(self, kmer):
        # hashes = []
        for i in range(self.hashFunctionCount):
            hashValue = mmh3.hash(kmer, i) % self.BFsize
            # hashes.append(hashValue)
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
