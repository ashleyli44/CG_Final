"""
Make the tree
"""
from queue import Queue
from node import Node


class Tree:

    def __init__(self, query_threshold):
        self.root = None
        self.query_threshold = query_threshold

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

        parent_bloom_filter = bloom_factory.new_instance()  # TODO: integrate with bloom factory
        parent_bloom_filter.union(curr_node.bloom_filter)
        parent_bloom_filter.union(new_node.bloom_filter)

        new_left_node = Node(curr_node.bloom_filter)
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

            if curr_match_count / float(len(kmers) + 1) >= self.query_threshold:  # +1 to avoid div by 0
                if curr_node.is_leaf():
                    potential_matches.append(curr_node)
                else:
                    bfs_queue.put(curr_node.left)
                    bfs_queue.put(curr_node.right)
