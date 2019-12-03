"""
Bloom filter factory
"""

from bitarray import bitarray
from bloom_filter import BloomFilter

class BloomFilterFactory:
    @staticmethod
    def new_instance(name=""):
        m = 1000
        h = 1
        return BloomFilter(name, bitarray(m), h)
