"""Union-find (disjoint-set) primitive shared by the clustering scripts."""


class UnionFind:
    """Disjoint-set with path compression over arbitrary hashable elements."""

    def __init__(self, elements):
        self.parent = {e: e for e in elements}

    def find(self, x):
        parent = self.parent
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb

    def groups(self):
        """Return ``{root: [members...]}`` in insertion order of members."""
        out = {}
        for e in self.parent:
            out.setdefault(self.find(e), []).append(e)
        return out