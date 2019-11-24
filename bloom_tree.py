"""
Make the tree
"""
from queue import Queue

from bloom_factory import BloomFilterFactory
from bloom_node import BloomNode


class BloomTree:

    def __init__(self, query_threshold):
        self.root = None
        self.query_threshold = query_threshold  # 0 < threshold < 1

    def add(self, new_node):
        if self.root is None:
            self.root = new_node
            return

        curr_node = self.root
        while not curr_node.is_leaf():
            curr_node.bloom_filter.union(new_node.bloom_filter)

            if curr_node.left is None:
                curr_node.left = new_node
                return

            if curr_node.right is None:
                curr_node.right = new_node
                return

            leftDistance = curr_node.left.bloom_filter.distanceTo(new_node.bloom_filter)
            rightDistance = curr_node.right.bloom_filter.distanceTo(new_node.bloom_filter)

            if rightDistance < leftDistance:
                curr_node = curr_node.right
            else:
                curr_node = curr_node.left

        parent_bloom_filter = BloomFilterFactory.new_instance()
        parent_bloom_filter.union(curr_node.bloom_filter)
        parent_bloom_filter.union(new_node.bloom_filter)

        new_left_node = BloomNode(curr_node.bloom_filter)
        curr_node.bloom_filter = parent_bloom_filter

        curr_node.left = new_left_node
        curr_node.right = new_node

    def query(self, kmers):
        if self.root is None:
            return None

        potential_matches = []
        bfs_queue = Queue()
        bfs_queue.put(self.root)

        while not bfs_queue.empty():
            curr_node = bfs_queue.get()
            if curr_node is None:
                continue

            curr_match_count = 0

            for kmer in kmers:
                if curr_node.bloom_filter.contains(kmer):
                    curr_match_count += 1

            if float(curr_match_count + 1) / float(len(kmers) + 2) >= self.query_threshold:  # +1/2 to avoid div by 0
                if curr_node.is_leaf():
                    potential_matches.append(curr_node)
                else:
                    bfs_queue.put(curr_node.left)
                    bfs_queue.put(curr_node.right)
