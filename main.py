from bloom_factory import BloomFilterFactory
from bloom_filter import BloomFilter
from bloom_node import BloomNode
from bloom_tree import BloomTree
from bitarray import bitarray
def main():
    tree = BloomTree(0.5);

    arr1 = bitarray(100)
    arr1.setall(False)
    bf1 = BloomFilter("a", arr1, 4)
    bf1.add("hello")
    arr2 = bitarray(100)
    arr2.setall(False)
    bf2 = BloomFilter("b", arr2, 4)
    bf2.add("potato")
    arr3 = bitarray(100)
    arr3.setall(False)
    bf3 = BloomFilter("c", arr3, 4)
    bf3.add("tomato")
    arr4 = bitarray(100)
    arr4.setall(False)
    bf4 = BloomFilter("d", arr4, 4)
    bf4.add("blah")
    """
    first = BloomNode(BloomFilter("a", tru_arr, 5))
    second = BloomNode(BloomFilter("b", fal_arr, 5))
    third = BloomNode(BloomFilter("c", tru_arr, 5))
    fourth = BloomNode(BloomFilter("d", fal_arr, 5))
    """
    first = BloomNode(bf1)
    second = BloomNode(bf2)
    third = BloomNode(bf3)
    fourth = BloomNode(bf4)
    tree.add(first)
    tree.add(second)
    tree.add(third)
    tree.add(fourth)

    tree_str = tree.treeToString(tree.root)
    print(tree_str)
    kmer_list = tree.query(["blah", "hello"])
    for i in kmer_list:
        print(i.bloom_filter.name)

main()
