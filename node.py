"""
TreeNode that holds a Bloom Filter
"""


class Node:

    def __init__(self, bloom_filter):
        self.bloom_filter = bloom_filter
        self.left = None
        self.right = None

    def is_leaf(self):
        return self.left is None and self.right is None
