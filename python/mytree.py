class LocalIndexer:
    def __init__(self, max_id):
        self.max_id = max_id
        self.stamp = [0] * (max_id + 1)
        self.localIdx = [0] * (max_id + 1)
        self.cur_stamp = 0

    def start_leaf(self, leaf_nodes):
        """Initialize a new leaf context with given global node IDs."""
        self.cur_stamp += 1
        for idx, node in enumerate(leaf_nodes):
            self.stamp[node] = self.cur_stamp
            self.localIdx[node] = idx

    def local_index(self, node):
        """Return the local index of `node` in the current leaf, or None if absent."""
        if 0 <= node <= self.max_id and self.stamp[node] == self.cur_stamp:
            return self.localIdx[node]
        return None

# Demonstration
indexer = LocalIndexer(max_id=10)

leaves = [
    [1, 5, 7],
    [2, 5, 9],
    [1, 2, 3, 9,10]
]

for i, leaf in enumerate(leaves, 1):
    indexer.start_leaf(leaf)


for i, leaf in enumerate(leaves, 1):
    print(f"Leaf {i}: {leaf}")
    for node in range(0, indexer.max_id + 1):
        local = indexer.local_index(node)
        if local is not None:
            print(f"  Global {node} -> Local {local}")
    print()

print("Final local indices:", indexer.localIdx)
print("Final stamps:", indexer.stamp)
print("Current stamp:", indexer.cur_stamp)
