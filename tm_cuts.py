import itertools, random, time
from my_enum.tm_library import get_bit, invert_truth_table

# Enumerates K-feasible cuts for every AIG node as the mapper's first step
# Each cut stores its root, leaves, absorbed inside nodes, and truth table for library matching and replication tracking.
# Builds cuts by merging fanin cuts, filtering invalid or redundant combinations, and optionally limiting the number retained per node.
# Also has cut validation

# AIG literals: literal = 2 * node + 1 if the edge is inverted
def lit_node(literal):
    return literal // 2

def lit_inv(literal):
    return literal % 2 == 1

class Cut:
    def __init__(self, root, leaves, inside, truth_table):
        self.root = root                   # the node this cut belongs to
        self.leaves = leaves               # tuple of node ids
        self.leafset = frozenset(leaves)   # the same leaves as a set
        self.inside = inside
        self.truth_table = truth_table     # root function of the leaves

    def __repr__(self):
        return f"Cut(root={self.root}, leaves={self.leaves}, inside={sorted(self.inside)}, truth_table={self.truth_table:#x})"

# ---------------------------------------------------------------------------- enumeration

def enumerate_cuts(aig, K=4, max_cuts=None, keep_shared_leaves=False):
    # produces cuts

    stats = dict(pairs=0, too_many_leaves=0, read_and_absorbed=0, repeated=0, redundant=0,
                 over_limit=0, limited_nodes=0)
    start = time.perf_counter()
    cuts = []
    for node in range(aig.n_nodes()): # node ids are in topological order
        trivial_cut = Cut(node, (node,), frozenset(), leaf_truth_table(0, 1))
        if aig.kind[node] != "and":
            cuts.append([trivial_cut]) # inputs and the constant: only the trivial cut
            continue
        fanin_a, fanin_b = lit_node(aig.fanin0[node]), lit_node(aig.fanin1[node])

        # merge every cut of fanin_a with every cut of fanin_b
        candidates = {}
        for cut_a in cuts[fanin_a]:
            for cut_b in cuts[fanin_b]:
                stats["pairs"] += 1
                leaf_set = cut_a.leafset | cut_b.leafset
                if len(leaf_set) > K:
                    stats["too_many_leaves"] += 1 # not K-feasible
                    continue
                inside_set = cut_a.inside | cut_b.inside | {node}
                if not leaf_set.isdisjoint(inside_set):
                    stats["read_and_absorbed"] += 1
                    continue
                if leaf_set in candidates:
                    stats["repeated"] += 1
                    continue
                candidates[leaf_set] = inside_set

        # remove redundant cuts, smallest first
        kept = []
        for leaf_set in sorted(candidates, key=len):
            if any(makes_redundant(aig, other, candidates[other], leaf_set, keep_shared_leaves)
                   for other in kept):
                stats["redundant"] += 1
            else:
                kept.append(leaf_set)

        # fewest leaves first
        kept.sort(key=lambda leaf_set: (len(leaf_set), sorted(leaf_set)))
        if max_cuts is not None and len(kept) > max_cuts:
            stats["over_limit"] += len(kept) - max_cuts
            stats["limited_nodes"] += 1
            kept = kept[:max_cuts]

        # build the Cut objects
        node_cuts = [trivial_cut]
        for leaf_set in kept:
            leaves = tuple(sorted(leaf_set))
            inside = candidates[leaf_set]
            node_cuts.append(Cut(node, leaves, inside, cut_function(aig, node, leaves, inside)))
        cuts.append(node_cuts)
    stats["seconds"] = time.perf_counter() - start
    return cuts, stats

def makes_redundant(aig, smaller_leaves, smaller_inside, leaves, keep_shared_leaves=False):
    # True if the cut with leaves smaller_leaves makes the cut with leaves `leaves` redundant
    if not smaller_leaves < leaves: # < on sets means "proper subset"
        return False
    if not keep_shared_leaves:
        return True

    # the smaller cut computes the extra leaves inside itself; that is only free when
    # nothing outside the smaller cut reads them
    for extra_leaf in leaves - smaller_leaves:
        if aig.po_refs[extra_leaf] > 0:
            return False # drives an output, so it keeps its own gate
        for reader in aig.fanouts[extra_leaf]:
            if reader not in smaller_inside:
                return False # read outside: absorbing it would replicate it
    return True

def leaf_truth_table(leaf, n_leaves):
    table = 0
    for row in range(2 ** n_leaves):
        if get_bit(row, leaf) == 1:
            table += 2 ** row
    return table

def cut_function(aig, root, leaves, inside):
    n_leaves = len(leaves)
    table = {0: 0}
    for i, leaf in enumerate(leaves):
        table[leaf] = leaf_truth_table(i, n_leaves)
    for node in sorted(inside): # sorted ids = topological order
        literal_a, literal_b = aig.fanin0[node], aig.fanin1[node]
        table_a = table[lit_node(literal_a)]
        if lit_inv(literal_a): # inverted edge: flip every row
            table_a = invert_truth_table(table_a, n_leaves)
        table_b = table[lit_node(literal_b)]
        if lit_inv(literal_b):
            table_b = invert_truth_table(table_b, n_leaves)
        table[node] = and_truth_tables(table_a, table_b, n_leaves)
    return table[root]

def and_truth_tables(table_a, table_b, n_leaves):
    result = 0
    for row in range(2 ** n_leaves):
        if get_bit(table_a, row) == 1 and get_bit(table_b, row) == 1:
            result += 2 ** row
    return result

# ---------------------------------------------------------------------------- checks

def check_cuts(aig, cuts, K, keep_shared_leaves=False, rounds=2, width=256, seed=565):
    # verify every cut

    problems = []

    def name(node):
        return aig.names[node]

    def describe(cut): # for example "u cut {a b c}"
        return f"{name(cut.root)} cut {{{' '.join(name(leaf) for leaf in cut.leaves)}}}"

    for node in range(aig.n_nodes()):
        node_cuts = cuts[node]

        # every node start with its trivial cut {node}; inputs and the constant have nothing else
        if not node_cuts or node_cuts[0].leaves != (node,) or node_cuts[0].inside:
            problems.append(f"{name(node)}: the first cut is not the trivial cut")
            continue
        if aig.kind[node] != "and" and len(node_cuts) > 1:
            problems.append(f"{name(node)}: an input or constant has non-trivial cuts")

        # 1. structure of each cut
        for cut in node_cuts[1:]:
            if cut.root != node or node not in cut.inside:
                problems.append(f"{describe(cut)}: its root is not among its inside nodes")
            if not 1 <= len(cut.leafset) <= K or len(cut.leafset) != len(cut.leaves):
                problems.append(f"{describe(cut)}: has {len(cut.leaves)} leaves, K = {K}")
            if cut.leafset & cut.inside:
                problems.append(f"{describe(cut)}: a node is both a leaf and inside")
            if any(aig.kind[inside_node] != "and" for inside_node in cut.inside):
                problems.append(f"{describe(cut)}: an input is listed as inside")
                continue

            # what the inside nodes read from outside the cut must be exactly the leavs
            read_from_outside = set()
            for inside_node in cut.inside:
                for literal in (aig.fanin0[inside_node], aig.fanin1[inside_node]):
                    if lit_node(literal) not in cut.inside:
                        read_from_outside.add(lit_node(literal))
            if read_from_outside != cut.leafset:
                read_names = " ".join(sorted(name(x) for x in read_from_outside))
                problems.append(f"{describe(cut)}: its inside nodes read {{{read_names}}}")

            # every inside node except the root must be read by another inside node
            for inside_node in cut.inside - {node}:
                if not any(reader in cut.inside for reader in aig.fanouts[inside_node]):
                    problems.append(f"{describe(cut)}: an inside node does not feed the root")
                    break

        # 2. no repeated or redundant cuts
        for position, cut in enumerate(node_cuts[1:]):
            for other_position, other in enumerate(node_cuts[1:]):
                if other_position == position:
                    continue                        # compare by position, so a cut listed twice is caught
                if other.leafset == cut.leafset or makes_redundant(
                        aig, other.leafset, other.inside, cut.leafset, keep_shared_leaves):
                    problems.append(f"{describe(cut)}: repeated, or made redundant by {describe(other)}")
                    break

    # 3. truth tables, checked on random input patterns (one pattern per bit of a `width`-bit number)
    random_generator = random.Random(seed)
    all_ones = 2 ** width - 1
    for _ in range(rounds):
        values = simulate_nodes(aig, random_generator, width)
        for node in aig.ands():
            for cut in cuts[node][1:]:
                leaf_values = [values[leaf] for leaf in cut.leaves]
                if apply_truth_table(cut.truth_table, leaf_values, all_ones) != values[node]:
                    problems.append(f"{describe(cut)}: truth table {cut.truth_table:#x} gives the wrong value")
    return problems

def simulate_nodes(aig, random_generator, width):
    all_ones = 2 ** width - 1 # 1 in every pattern
    values = [0] * aig.n_nodes()
    for input_node in aig.pis:
        values[input_node] = random_generator.getrandbits(width)
    for node in aig.ands():
        literal_a, literal_b = aig.fanin0[node], aig.fanin1[node]
        value_a = values[lit_node(literal_a)]
        if lit_inv(literal_a):
            value_a = all_ones - value_a
        value_b = values[lit_node(literal_b)]
        if lit_inv(literal_b):
            value_b = all_ones - value_b
        values[node] = value_a & value_b
    return values

def apply_truth_table(truth_table, leaf_values, all_ones):
    output = 0
    for row in range(2 ** len(leaf_values)):
        if get_bit(truth_table, row) == 0:
            continue # the output is 0 in this row
        matching = all_ones # patterns that match this row so far: all of them
        for leaf, value in enumerate(leaf_values):
            if get_bit(row, leaf) == 1:
                matching = matching & value # kep the patterns where this leaf is 1
            else:
                matching = matching & (all_ones - value)
        output = output | matching
    return output


def brute_force_cuts(aig, K, keep_shared_leaves=False):
    result = {}
    for node in aig.ands():
        nodes_below = sorted(fanin_cone(aig, node))
        real_cuts = {} # leaf set -> inside set
        for size in range(1, K + 1):
            for leaves in itertools.combinations(nodes_below, size):   # every set of `size` nodes
                leaf_set = frozenset(leaves)
                inside = cone(aig, node, leaf_set)
                if inside is not None: # not a real cut
                    real_cuts[leaf_set] = inside

        # keep only the cuts that no other real cut makes redundant
        result[node] = {}
        for leaf_set, inside in real_cuts.items():
            redundant = False
            for other in real_cuts:
                if makes_redundant(aig, other, real_cuts[other], leaf_set, keep_shared_leaves):
                    redundant = True
                    break
            if not redundant:
                result[node][leaf_set] = inside
    return result

def fanin_cone(aig, node):
    below = set()
    to_visit = [node] # nodes whose inputs still have to be looked at
    while to_visit:
        current = to_visit.pop()
        if aig.kind[current] == "and": # inputs and the constant have nothing below
            for fanin in (lit_node(aig.fanin0[current]), lit_node(aig.fanin1[current])):
                if fanin not in below:
                    below.add(fanin)
                    to_visit.append(fanin)
    return below

def cone(aig, node, leaf_set):
    inside = set()
    used_leaves = set()
    to_visit = [node]
    while to_visit:
        current = to_visit.pop()
        if current in leaf_set:
            used_leaves.add(current)
        elif current not in inside:
            if aig.kind[current] != "and":
                return None
            inside.add(current)
            to_visit += [lit_node(aig.fanin0[current]), lit_node(aig.fanin1[current])]
    if used_leaves != leaf_set:
        return None # unused
    return frozenset(inside)

# ---------------------------------------------------------------------------- output

def print_summary(aig, cuts, stats, K):
    and_nodes = aig.ands()
    cut_counts = [len(cuts[node]) - 1 for node in and_nodes]
    total = sum(cut_counts)
    by_size = {}
    for node in and_nodes:
        for cut in cuts[node][1:]:
            by_size[len(cut.leaves)] = by_size.get(len(cut.leaves), 0) + 1

    circuit = aig.name.replace("\\", "/").split("/")[-1] # drop any folder path from the name
    print(f"{circuit}: {len(aig.pis)} inputs, {len(and_nodes)} AND nodes, K = {K}")
    if and_nodes:
        print(f"cuts      : {total} ({total / len(and_nodes):.1f} per AND node, at most {max(cut_counts)})")
    print("by leaves : " + ", ".join(f"{size} leaves: {by_size[size]}" for size in sorted(by_size)))
    print(f"pairs     : {stats['pairs']} tried; skipped {stats['too_many_leaves']} with too many leaves, "
          f"{stats['read_and_absorbed']} read-and-absorbed, {stats['repeated']} repeated; "
          f"removed {stats['redundant']} redundant")
    if stats["limited_nodes"]:
        print(f"limit     : {stats['over_limit']} cuts dropped at {stats['limited_nodes']} nodes")
    print(f"time      : {stats['seconds']:.2f} s")


def print_cuts(aig, cuts, name):
    node = aig.names.index(name)
    print(f"cuts of {name} (besides its trivial cut {{{name}}}):")
    for number, cut in enumerate(cuts[node][1:], start=1):
        leaves = " ".join(aig.names[leaf] for leaf in cut.leaves)
        inside = " ".join(aig.names[inside_node] for inside_node in sorted(cut.inside))
        rows = format(cut.truth_table, f"0{2 ** len(cut.leaves)}b")
        print(f"  {number:2d}. leaves {{{leaves}}}   inside {{{inside}}}   truth table {rows}")