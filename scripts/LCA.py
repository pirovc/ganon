"""LCA.py

Range minimization and tree least common ancestor data structures
with linear space and preprocessing time, and constant query time,
from Bender and Farach-Colton, "The LCA Problem Revisited",
Proc. LATIN 2000 (pp.88-94), http://www.cs.sunysb.edu/~bender/pub/lca.ps

Some experimentation would be needed to determine how large a query
range needs to be to make this faster than computing the min of the range
directly, and how much input data is needed to make the linear space
version pay off compared to the much simpler LogarithmicRangeMin that
it uses as a subroutine.

D. Eppstein, November 2003.
"""

import unittest,random
from collections import defaultdict
#from UnionFind import UnionFind

"""UnionFind.py

Union-find data structure. Based on Josiah Carlson's code,
http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/215912
with significant additional changes by D. Eppstein.
"""
class UnionFind:
    """Union-find data structure.

    Each unionFind instance X maintains a family of disjoint sets of
    hashable objects, supporting the following two methods:

    - X[item] returns a name for the set containing the given item.
      Each set is named by an arbitrarily-chosen one of its members; as
      long as the set remains unchanged it will keep the same name. If
      the item is not yet part of a set in X, a new singleton set is
      created for it.

    - X.union(item1, item2, ...) merges the sets containing each item
      into a single larger set.  If any item is not yet part of a set
      in X, it is added to X as one of the members of the merged set.
    """

    def __init__(self):
        """Create a new empty union-find structure."""
        self.weights = {}
        self.parents = {}

    def __getitem__(self, object):
        """Find and return the name of the set containing the object."""

        # check for previously unknown object
        if object not in self.parents:
            self.parents[object] = object
            self.weights[object] = 1
            return object

        # find path of objects leading to the root
        path = [object]
        root = self.parents[object]
        while root != path[-1]:
            path.append(root)
            root = self.parents[root]

        # compress the path and return
        for ancestor in path:
            self.parents[ancestor] = root
        return root
        
    def __iter__(self):
        """Iterate through all items ever found or unioned by this structure."""
        return iter(self.parents)

    def union(self, *objects):
        """Find the sets containing the objects and merge them all."""
        roots = [self[x] for x in objects]
        heaviest = max([(self.weights[r],r) for r in roots])[1]
        for r in roots:
            if r != heaviest:
                self.weights[heaviest] += self.weights[r]
                self.parents[r] = heaviest
				
# 2to3 compatibility
try:
    xrange
except:
    xrange = range

def _decodeSlice(self,it):
    """Work around removal of __getslice__ in Python 3"""
    if type(it) != slice:
        raise ValueError("Can only access LCA object by slice notation")
    left,right,stride = it.indices(len(self))
    if stride != 1:
        raise ValueError("Stride not permitted in LCA")
    return left,right

class RangeMin:
    """If X is any list, RangeMin(X)[i:j] == min(X[i:j]).
    Initializing RangeMin(X) takes time and space linear in len(X),
    and querying the minimum of a range takes constant time per query.
    """

    def __init__(self,X):
        """Set up structure with sequence X as data.
        Uses an LCA structure on a Cartesian tree for the input."""
        self._data = list(X)
        if len(X) > 1:
            big = list(map(max, self._ansv(False), self._ansv(True)))
            parents = {i:big[i][1] for i in range(len(X)) if big[i]}
            self._lca = LCA(parents)

    def __getitem__(self,it):
        """When called by X[left:right], return min(X[left:right])."""
        left,right = _decodeSlice(self,it)
        if right <= left:
            return None     # empty range has no minimum
        return self._data[self._lca(left,right-1)]

    def __len__(self):
        """How much data do we have?  Needed for negative index in slice."""
        return len(self._data)

    def _ansv(self,reversed):
        """All nearest smaller values.
        For each x in the data, find the value smaller than x in the closest
        position to the left of x (if not reversed) or to the right of x
        (if reversed), and return list of pairs (smaller value,position).
        Due to our use of positions as a tie-breaker, values equal to x
        count as smaller on the left and larger on the right.
        """
        stack = [(min(self._data),-1)]   # protect stack top with sentinel
        output = [0]*len(self._data)
        for xi in _pairs(self._data,reversed):
            while stack[-1] > xi:
                stack.pop()
            output[xi[1]] = stack[-1]
            stack.append(xi)
        return output

    def _lca(self,first,last):
        """Function to replace LCA when we have too little data."""
        return 0

class RestrictedRangeMin:
    """Linear-space RangeMin for integer data obeying the constraint
        abs(X[i]-X[i-1])==1.
    We don't actually check this constraint, but results may be incorrect
    if it is violated.  For the use of this data structure from LCA, the
    data are actually pairs rather than integers, but the minima of
    all ranges are in the same positions as the minima of the integers
    in the first positions of each pair, so the data structure still works.
    """
    def __init__(self,X):
        # Compute parameters for partition into blocks.
        # Position i in X becomes transformed into
        # position i&self._blockmask in block i>>self.blocklen
        self._blocksize = _log2(len(X))//2
        self._blockmask = (1 << self._blocksize) - 1
        blocklen = 1 << self._blocksize

        # Do partition into blocks, find minima within
        # each block, prefix minima in each block,
        # and suffix minima in each block
        blocks = []             # map block to block id
        ids = {}                # map block id to PrecomputedRangeMin
        blockmin = []           # map block to min value
        self._prefix = [None]   # map data index to prefix min of block
        self._suffix = []       # map data index to suffix min of block
        for i in range(0,len(X),blocklen):
            XX = X[i:i+blocklen]
            blockmin.append(min(XX))
            self._prefix += PrefixMinima(XX)
            self._suffix += PrefixMinima(XX,reversed=True)
            blockid = len(XX) < blocklen and -1 or self._blockid(XX)
            blocks.append(blockid)
            if blockid not in ids:
                ids[blockid] = PrecomputedRangeMin(_pairs(XX))
        self._blocks = [ids[b] for b in blocks]

        # Build data structure for interblock queries
        self._blockrange = LogarithmicRangeMin(blockmin)
        self._data = list(X)

    def __len__(self):
        """How much data do we have?  Needed for negative index in slice."""
        return len(self._data)

    def __getitem__(self,it):
        """When called by X[left:right], return min(X[left:right])."""
        left,right = _decodeSlice(self,it)
        firstblock = left >> self._blocksize
        lastblock = (right - 1) >> self._blocksize
        if firstblock == lastblock:
            i = left & self._blockmask
            position = self._blocks[firstblock][i:i+right-left][1]
            return self._data[position + (firstblock << self._blocksize)]
        else:
            best = min(self._suffix[left], self._prefix[right])
            if lastblock > firstblock + 1:
                best = min(best, self._blockrange[firstblock+1:lastblock])
            return best

    def _blockid(self,XX):
        """Return value such that all blocks with the same
        pattern of increments and decrements get the same id.
        """
        blockid = 0
        for i in range(1,len(XX)):
            blockid = blockid*2 + (XX[i] > XX[i-1])
        return blockid

class PrecomputedRangeMin:
    """RangeMin solved in quadratic space by precomputing all solutions."""

    def __init__(self,X):
        self._minima = [PrefixMinima(X[i:]) for i in range(len(X))]

    def __getitem__(self,it):
        """When called by X[left:right], return min(X[left:right])."""
        left,right = _decodeSlice(self,it)
        return self._minima[left][right-left-1]

    def __len__(self):
        return len(self._minima)

class LogarithmicRangeMin:
    """RangeMin in O(n log n) space and constant query time."""

    def __init__(self,X):
        """Compute min(X[i:i+2**j]) for each possible i,j."""
        self._minima = m = [list(X)]
        for j in range(_log2(len(X))):
            m.append(list(map(min, m[-1][:-1<<j], m[-1][1<<j:])))

    def __getitem__(self,it):
        """When called by X[left:right], return min(X[left:right])."""
        left,right = _decodeSlice(self,it)
        j = _logtable[right-left]
        row = self._minima[j]
        return min(row[left],row[right-2**j])

    def __len__(self):
        return len(self._minima[0])

class LCA:
    """Structure for finding least common ancestors in trees.
    Tree nodes may be any hashable objects; a tree is specified
    by a dictionary mapping nodes to their parents.
    LCA(T)(x,y) finds the LCA of nodes x and y in tree T.
    """
    def __init__(self, parent, RangeMinFactory = RestrictedRangeMin):
        """Construct LCA structure from tree parent relation."""
        children = defaultdict(list)
        for x in parent:
            children[parent[x]].append(x)
        root = [x for x in children if x not in parent]
        if len(root) != 1:
            raise ValueError("LCA input is not a tree")
        
        levels = []
        self._representatives = {}
        self._visit(children,levels,root[0],0)
        if [x for x in parent if x not in self._representatives]:
            raise ValueError("LCA input is not a tree")
        self._rangemin = RangeMinFactory(levels)

    def __call__(self,*nodes):
        """Find least common ancestor of a set of nodes."""
        r = [self._representatives[x] for x in nodes]
        return self._rangemin[min(r):max(r)+1][1]

    def _visit(self,children,levels,node,level):
        """Perform Euler traversal of tree."""
        self._representatives[node] = len(levels)
        pair = (level,node)
        levels.append(pair)
        for child in children[node]:
            self._visit(children,levels,child,level+1)
            levels.append(pair)

class OfflineLCA(defaultdict):
    """Find LCAs of all pairs in a given sequence, using Union-Find."""

    def __init__(self,parent,pairs):
        """Set up to find LCAs of pairs in tree defined by parent.
        LCA of any given pair x,y can then be found by self[x][y].
        However unlike the online LCA structure we can not find LCAs
        of pairs that are not supplied to us at startup time.
        """

        # set up dictionary where answers get stored
        defaultdict.__init__(self,dict)
        for u,v in pairs:
            self[u][v] = self[v][u] = None

        # set up data structure for finding node ancestor on search path
        # self.descendants forms a collection of disjoint sets,
        #    one set for the descendants of each search path node.
        # self.ancestors maps disjoint set ids to the ancestors themselves.
        self.descendants = UnionFind()
        self.ancestors = {}

        # invert the parent relationship so we can traverse the tree
        self.children = defaultdict(list)
        for x,px in parent.items():
            self.children[px].append(x)
        root = [x for x in self.children if x not in parent]
        if len(root) != 1:
            raise ValueError("LCA input is not a tree")

        # initiate depth first traversal
        self.visited = set()
        self.traverse(root[0])

    def traverse(self,node):
        """Perform depth first traversal of tree."""
        self.ancestors[self.descendants[node]] = node
        for child in self.children[node]:
            self.traverse(child)
            self.descendants.union(child,node)
            self.ancestors[self.descendants[node]] = node
        self.visited.add(node)
        for query in self[node]:
            if query in self.visited:
                lca = self.ancestors[self.descendants[query]]
                self[node][query] = self[query][node] = lca

# Various utility functions

def PrefixMinima(X,reversed=False):
    """Compute table of prefix minima
    (or suffix minima, if reversed=True) of list X.
    """
    current = None
    output = [None]*len(X)
    for x,i in _pairs(X,reversed):
        if current is None:
            current = x
        else:
            current = min(current,x)
        output[i] = current
    return output

def _pairs(X,reversed=False):
    """Return pairs (x,i) for x in list X, where i is
    the index of x in the data, in forward or reverse order.
    """
    if reversed:
        indices = range(len(X)-1,-1,-1)
    else:
        indices = range(len(X))
    return [(X[i],i) for i in indices]

_logtable = [None,0]
def _log2(n):
    """Make table of logs reach up to n and return floor(log_2(n))."""
    while len(_logtable) <= n:
        _logtable.extend([1+_logtable[-1]]*len(_logtable))
    return _logtable[n]

# if run as "python LCA.py", run tests on random data
# and check that RangeMin's results are correct.

class RandomRangeMinTest(unittest.TestCase):
    def testRangeMin(self):
        for trial in range(20):
            data = [random.choice(xrange(1000000))
                    for i in range(random.randint(1,100))]
            R = RangeMin(data)
            for sample in range(100):
                i = random.randint(0,len(data)-1)
                j = random.randint(i+1,len(data))
                self.assertEqual(R[i:j],min(data[i:j]))

class LCATest(unittest.TestCase):
    parent = {'b':'a','c':'a','d':'a','e':'b','f':'b','g':'f','h':'g','i':'g'}
    lcas = {
        ('a','b'):'a',
        ('b','c'):'a',
        ('c','d'):'a',
        ('d','e'):'a',
        ('e','f'):'b',
        ('e','g'):'b',
        ('e','h'):'b',
        ('c','i'):'a',
        ('a','i'):'a',
        ('f','i'):'f',
    }

    def testLCA(self):
        L = LCA(self.parent)
        for k,v in self.lcas.items():
            self.assertEqual(L(*k),v)

    def testLogLCA(self):
        L = LCA(self.parent, LogarithmicRangeMin)
        for k,v in self.lcas.items():
            self.assertEqual(L(*k),v)

    def testOfflineLCA(self):
        L = OfflineLCA(self.parent, self.lcas.keys())
        for (p,q),v in self.lcas.items():
            self.assertEqual(L[p][q],v)

	
if __name__ == "__main__":
    unittest.main()   


	
