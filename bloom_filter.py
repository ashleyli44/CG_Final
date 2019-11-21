'''
Bloom filter class
'''
import bitarray


class BoomFilter:
    def __init__(self, str_id, b_vect):
        self.name = str_id #string name
        self.bvector = b_vect #bitarray vector that represents the bloom filter
    
    def compareTo(self, bf):
        '''
            Calculates the hamming distance b/w itself and another bf
            @param bf is the bloom filter object to pass in --> assume the length is the same
            @return is the integer representing the hamming distance
            w/ 0 being the minimum and len(bf) being the maximum
        '''
        #exclusive OR and then calculates the number of set bits
        return (self.bvector ^ bf.vector).count(True)
                 
