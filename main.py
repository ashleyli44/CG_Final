from bloom_factory import BloomFilterFactory
from bloom_filter import BloomFilter
from bloom_node import BloomNode
from bloom_tree import BloomTree
from bitarray import bitarray
def main():
    tree = BloomTree(0.5);
    first = BloomNode(BloomFilter("a", bitarray(1000), 5))
    second = BloomNode(BloomFilter("b", bitarray(1000), 5))
    third = BloomNode(BloomFilter("c", bitarray(1000), 5))
    fourth = BloomNode(BloomFilter("d", bitarray(1000), 5))
    tree.add(first)
    tree.add(second)
    tree.add(third)
    tree.add(fourth)

main()

print("hello2")