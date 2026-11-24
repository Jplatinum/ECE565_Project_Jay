import itertools

"""
Reads the genlib cell library and builds a truth-table lookup for matching AIG cuts to cells.
Stores cell area, delay parameters, and all possible pin orderings and input/output inversion choices.
Provides truth-table utilities and functions for parsing genlib Boolean expressions.
Main interface: read_genlib() loads the library, and Library.match_table() retrieves matching cells for a cut.
"""

# cell class that hold intrinsic details to the cell
# Pins will keep a list of the pin input names Example: ["a", "b"]
# truth table will map the matching for the checking simulation
class Cell:
    def __init__(self, name, area, pins, C, d, R, truth_table=0):
        self.name = name
        self.area = area
        self.pins = pins
        self.C = C
        self.d = d
        self.R = R
        self.truth_table = truth_table

    # to get number of fanins will be number of pins
    def n_inputs(self):
        return len(self.pins)

# track a set of cells which will be our library
class Library:
    def __init__(self, cells):
        self.cells = cells # list of Cell, in the order of the genlib file

    def cell(self, name):
        return next(c for c in self.cells if c.name == name)

    def max_inputs(self):
        return max(c.n_inputs() for c in self.cells)

    # builds the lookup table that matches cuts to cells
    def match_table(self):
        table = {}
        for cell in self.cells:
            n_pins = cell.n_inputs()
            if n_pins == 0:
                continue

            # every order in which the pins can be connected
            for pin_to_leaf in itertools.permutations(range(n_pins)):
                # every choice of which pin is inverted
                for choice in range(2 ** n_pins):
                    inverted_pins = [get_bit(choice, pin) == 1 for pin in range(n_pins)]
                    function = function_seen_by_leaves(cell, pin_to_leaf, inverted_pins)

                    # record output
                    for output_inverted in (False, True):
                        if output_inverted:
                            key = (n_pins, invert_truth_table(function, n_pins))
                        else:
                            key = (n_pins, function)
                        if key not in table:
                            table[key] = []
                        table[key].append(Match(cell, pin_to_leaf, inverted_pins, output_inverted))
        return table


class Match:
    def __init__(self, cell, pin_to_leaf, inverted_pins, output_inverted):
        self.cell = cell
        self.pin_to_leaf = pin_to_leaf
        self.inverted_pins = inverted_pins
        self.output_inverted = output_inverted

def get_bit(number, position):
    # helper functino
    return (number // 2 ** position) % 2

def invert_truth_table(truth_table, n_inputs):
    all_rows_one = 2 ** (2 ** n_inputs) - 1 # the truth table that is 1 in every row
    return all_rows_one - truth_table

def function_seen_by_leaves(cell, pin_to_leaf, inverted_pins):
    # truth table for the cell based on the cuts leaves
    n_pins = cell.n_inputs()
    function = 0
    for leaf_row in range(2 ** n_pins):
        cell_row = 0
        for pin in range(n_pins):
            pin_value = get_bit(leaf_row, pin_to_leaf[pin])
            if inverted_pins[pin]:
                pin_value = 1 - pin_value
            cell_row += pin_value * 2 ** pin
        cell_output = get_bit(cell.truth_table, cell_row)
        function += cell_output * 2 ** leaf_row
    return function

# One-character symbols in a genlib function:
SYMBOLS = "!*&+|^()'"
def _tokenize(text):
    # Split a genlib function such as "!(a*b)" into tokens ['!', '(', 'a', '*', 'b', ')']
    tokens = []
    position = 0
    while position < len(text):
        char = text[position]
        if char.isspace():
            position += 1 # spaces only separate tokens
        elif char in SYMBOLS:
            tokens.append(char) # a symbol is a token by itself
            position += 1
        else:
            end = position
            while end < len(text) and not text[end].isspace() and text[end] not in SYMBOLS:
                end += 1
            tokens.append(text[position:end])
            position = end
    return tokens

def _parse_expr(text):
    # Parse a genlib function into a tree of nested tuples
    tokens = _tokenize(text)
    position = 0

    def peek():
        if position < len(tokens):
            return tokens[position]
        return None

    def take():
        nonlocal position # lets this function change position
        position += 1
        return tokens[position - 1]

    def parse_or():
        tree = parse_xor()
        while peek() in ("+", "|"):
            take()
            tree = ("or", tree, parse_xor())
        return tree

    def parse_xor():
        tree = parse_and()
        while peek() == "^":
            take()
            tree = ("xor", tree, parse_and())
        return tree

    def parse_and():
        tree = parse_not()
        while peek() is not None and peek() not in ("+", "|", "^", ")"):
            if peek() in ("*", "&"):
                take()
            tree = ("and", tree, parse_not())
        return tree

    def parse_not():
        if peek() == "!":
            take()
            return ("not", parse_not())
        tree = parse_atom()
        while peek() == "'":
            take()
            tree = ("not", tree)
        return tree

    def parse_atom():
        token = take()
        if token == "(":
            tree = parse_or() # inside parentheses, start at weakest level
            if take() != ")":
                raise ValueError(f"missing ')' in {text!r}")
            return tree
        if token in ("CONST0", "CONST1"):
            return ("const", int(token[-1]))
        return ("var", token)

    tree = parse_or()                           # the whole function, starting at the weakest level
    if position != len(tokens):
        raise ValueError(f"cannot parse {text!r}")   # some tokens were left over
    return tree

def _vars(tree, pin_names):
    # collects pin names
    kind = tree[0]
    if kind == "var":
        name = tree[1]
        if name not in pin_names:
            pin_names.append(name)
    elif kind in ("not", "and", "or", "xor"):
        for operand in tree[1:]: # the operands left to right
            _vars(operand, pin_names)
    return pin_names

def _eval(tree, pin_values):
    kind = tree[0]
    if kind == "var":
        return pin_values[tree[1]] # the value given for this pin
    if kind == "const":
        return tree[1] # 0 or 1
    if kind == "not":
        return 1 - _eval(tree[1], pin_values) # 0 -> 1, 1 -> 0

    # and / or / xor: work out both operands first
    left = _eval(tree[1], pin_values)
    right = _eval(tree[2], pin_values)
    if kind == "and":
        return 1 if left == 1 and right == 1 else 0
    if kind == "or":
        return 1 if left == 1 or right == 1 else 0
    return 1 if left != right else 0 # xor

def read_genlib(path):
    # Read a genlib file  into a Library in the form of our defined cells

    with open(path) as file:
        # drop comments (from # to the end of a line)
        text = "\n".join(line.split("#")[0] for line in file)
    words = text.split()

    cells = []
    position = 0 # index of the next word to read
    while position < len(words):
        if words[position] != "GATE":
            raise ValueError(f"{path}: expected GATE, got {words[position]!r}")
        name = words[position + 1]
        area = float(words[position + 2])
        position += 3

        # the function can contain spaces, collect words ends with ';'
        function_words = []
        while True:
            function_words.append(words[position])
            position += 1
            if function_words[-1].endswith(";"):
                break
        function_text = " ".join(function_words).rstrip(";")
        function_text = function_text.split("=", 1)[1].strip()
        tree = _parse_expr(function_text)
        pins = _vars(tree, []) # pin order = order of first appearance

        pin_numbers = {}
        while position < len(words) and words[position] == "PIN":
            pin_name = words[position + 1]
            input_load = float(words[position + 3])
            rise_block, rise_fanout = float(words[position + 5]), float(words[position + 6])
            fall_block, fall_fanout = float(words[position + 7]), float(words[position + 8])
            pin_numbers[pin_name] = (input_load, max(rise_block, fall_block), max(rise_fanout, fall_fanout))
            position += 9

        # one C, d and R for the cell
        C = d = R = 0.0
        for pin in pins:
            if pin in pin_numbers:
                pin_C, pin_d, pin_R = pin_numbers[pin]
            elif "*" in pin_numbers:
                pin_C, pin_d, pin_R = pin_numbers["*"]
            else:
                raise ValueError(f"{path}: gate {name} has no PIN line for {pin}")
            C, d, R = max(C, pin_C), max(d, pin_d), max(R, pin_R)

        truth_table = 0
        for row in range(2 ** len(pins)):
            pin_values = {pin: get_bit(row, k) for k, pin in enumerate(pins)}
            truth_table += _eval(tree, pin_values) * 2 ** row

        cells.append(Cell(name, area, pins, C, d, R, truth_table))
    return Library(cells)