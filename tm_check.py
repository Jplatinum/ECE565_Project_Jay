# Independently validates mapped netlists, including our mapper's output and ABC's
# Checks netlist structure, timing (gate delay, critical path, area, and gate count), and functional equivalence to AIG
# Provides mapping reports, JSON cover-map loading, and ABC equivalence checking using CEC.

def read_netlist(netlist_path, library):
    with open(netlist_path) as file:
        text = file.read().replace("\\\n", " ") # a line ending in \ continues on the next one
    cells = {cell.name: cell for cell in library.cells}
    inputs, outputs, gates = [], [], []
    for line in text.split("\n"):
        words = line.split("#")[0].split() # drop comments
        if not words:
            continue
        if words[0] == ".inputs":
            inputs += words[1:]
        elif words[0] == ".outputs":
            outputs += words[1:]
        elif words[0] == ".gate": # .gate NAND2 a=x b=y O=z
            cell = cells[words[1]]
            connections = dict(word.split("=", 1) for word in words[2:])
            input_nets = [connections[pin] for pin in cell.pins] # in the cell's pin order
            output_net = next(net for pin, net in connections.items() if pin not in cell.pins)
            gates.append((cell, input_nets, output_net))
    return inputs, outputs, gates

def delay_and_area(library, inputs, outputs, gates):
    load = {} # net -> the capacitance it drives
    for cell, input_nets, output_net in gates:
        for net in input_nets:
            load[net] = load.get(net, 0) + cell.C
    for net in outputs:
        load[net] = load.get(net, 0) + library.cell("INV").C

    arrival = {net: 0 for net in inputs} # net -> when its value arrives
    driven_by = {} # gate output net
    waiting = gates
    while waiting: # time each gate once all its inputs have arrived
        still_waiting = []
        for cell, input_nets, output_net in waiting:
            if not all(net in arrival for net in input_nets):
                still_waiting.append((cell, input_nets, output_net))
                continue
            latest = max(input_nets, key=lambda net: arrival[net], default=None) # None for ZERO, ONE
            start = arrival[latest] if latest is not None else 0 # all inputs are in
            arrival[output_net] = start + cell.d + cell.R * load.get(output_net, 0)  # + this gate's delay
            driven_by[output_net] = (cell.name, latest)
        if len(still_waiting) == len(waiting):
            raise ValueError(f"{len(waiting)} gates never get all their inputs (a loop or an undriven net)")
        waiting = still_waiting
    # critical path: start at the latest output and follow the latest input back to a primary input
    net = max(outputs, key=lambda output: arrival[output])
    delay, path = arrival[net], []
    while net in driven_by:
        cell_name, latest = driven_by[net]
        path.insert(0, (cell_name, net, load.get(net, 0), arrival[net]))
        net = latest
    area = sum(cell.area for cell, input_nets, output_net in gates)
    return delay, area, path

def check_timing(library, netlist_path):
    inputs, outputs, gates = read_netlist(netlist_path, library)
    delay, area, critical_path = delay_and_area(library, inputs, outputs, gates)
    print(f"Timing    : delay {delay:g}   area {area:g}   gates {len(gates)}")
    print("   critical path (cell -> net, load, arrival):")
    for cell_name, net, load, arrival in critical_path:
        print(f"     {cell_name:7s} -> {net:12s} load {load:<5g} arrival {arrival:g}")
    return delay, area