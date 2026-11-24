from tm_cuts import enumerate_cuts, check_cuts, lit_node, lit_inv

# Demonstrates cut enumeration on three small AIG circuits.
# Prints each AND node's cuts as leaf-name tuples and validates them using check_cuts.
# Includes examples of shared logic

class SmallAIG:
    def __init__(self, inputs, gates, outputs):
        self.names = ['const0'] + inputs # node 0 is the constant, then the inputs
        self.kind = ['const'] + ['pi'] * len(inputs)
        self.fanin0 = [0] * len(self.names)
        self.fanin1 = [0] * len(self.names)
        # for every gate append the information
        for gate, (signal_a, signal_b) in gates.items():
            self.fanin0.append(self.literal(signal_a))
            self.fanin1.append(self.literal(signal_b))
            self.names.append(gate)
            self.kind.append('and')
        self.pis = list(range(1, len(inputs) + 1))
        self.pos = [(output, self.literal(output)) for output in outputs]
        self.fanouts = [[] for _ in self.names] # the gates that read each node
        for node in self.ands():
            for fanin in (lit_node(self.fanin0[node]), lit_node(self.fanin1[node])):
                if node not in self.fanouts[fanin]:
                    self.fanouts[fanin].append(node)
        self.po_refs = [0] * len(self.names) # how many outputs each node drives
        for _, literal in self.pos:
            self.po_refs[lit_node(literal)] += 1

    def literal(self, signal):
        node = self.names.index(signal.lstrip('!'))
        return 2 * node + 1 if signal.startswith('!') else 2 * node

    def n_nodes(self):
        return len(self.names)

    def ands(self):
        return [node for node in range(self.n_nodes()) if self.kind[node] == 'and']

# ---------------------------------------------------------------------------- examples/runs

# Example 1: one shared node
def shared_node_example():
    return SmallAIG(['a', 'b', 'x', 'y'], {
        'u': ('a', 'b'), 'v': ('u', 'x'), 'w': ('u', 'y'),
    }, ['!v', '!w'])


# Example 2: the 8-node NAND2/INV tree from my HW1
def hw1_example():
    return SmallAIG(['x1', 'x2', 'x3', 'x4', 'x5', 'x6'], {
        'n7': ('x1', 'x2'), 'n8': ('x3', 'x4'), 'n5': ('n7', '!n8'),
        'n4': ('x5', 'x6'), 'n2': ('n5', '!n4'),
    }, ['n2'])


# Example 3: u, v and q are shared --> u's branches reconverge at y
def larger_example():
    return SmallAIG(['a', 'b', 'c', 'd', 'e', 'f'], {
        'u': ('a', 'b'), 'v': ('c', 'd'),
        'p': ('u', 'e'), 'q': ('u', 'f'),
        'r': ('p', 'v'), 's': ('q', '!v'),
        'y': ('r', 'q'), 'z': ('s', 'e'),
    }, ['y', 'z'])

# ---------------------------------------------------------------------------- demo output

def run_demo(name, aig, K=4):
    # Enumerate the cuts, show each gate with its cuts, then check every cut.
    cuts, stats = enumerate_cuts(aig, K)
    names = aig.names

    def signal(literal):
        return ('!' if lit_inv(literal) else '') + names[lit_node(literal)]

    print(f'\n{name}')
    print(f'Cut enumeration: K = {K} (largest cell in our library) + trivial cuts not shown')
    print('Shared nodes:', [names[node] for node in aig.ands() if len(aig.fanouts[node]) > 1])
    candidates = {}
    for node in aig.ands():
        candidates[names[node]] = [tuple(names[leaf] for leaf in cut.leaves) for cut in cuts[node][1:]]
        print(f'  {names[node]} = {signal(aig.fanin0[node])} AND {signal(aig.fanin1[node])}; '
              f'cuts: {candidates[names[node]]}')
    print(f"  Pairs tried: {stats['pairs']}, too many leaves: {stats['too_many_leaves']}, "
          f"read-and-absorbed: {stats['read_and_absorbed']}, redundant: {stats['redundant']}")
    problems = check_cuts(aig, cuts, K)
    print('  Check:', 'PASS, every cut verified' if not problems else problems)
    return candidates


if __name__ == '__main__':
    run_demo('Example 1: one shared node', shared_node_example())
    run_demo('Example 2: HW1 tree', hw1_example())
    run_demo('Example 3: several shared nodes and an inversion', larger_example())