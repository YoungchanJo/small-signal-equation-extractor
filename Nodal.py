from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction

from SolverKernel import AtomTable, ExprArena, SPolyArena, UnknownInfo, plain, nodeName


ground_names = {"0", "gnd"}
mos_rows = ("g", "d", "s", "b")
bjt_rows = ("c", "b", "e")
mos_cap_pairs = {"gs": ("g", "s"), "gd": ("g", "d"), "ds": ("d", "s")}
mos_body_cap_pairs = mos_cap_pairs | {"gb": ("g", "b"), "db": ("d", "b"), "sb": ("s", "b")}
bjt_cap_pairs = {"be": ("b", "e"), "bc": ("b", "c"), "ce": ("c", "e")}
unit_scale = {
    "f": Fraction(1, 10**15),
    "p": Fraction(1, 10**12),
    "n": Fraction(1, 10**9),
    "u": Fraction(1, 10**6),
    "m": Fraction(1, 10**3),
    "k": Fraction(10**3, 1),
    "meg": Fraction(10**6, 1),
    "g": Fraction(10**9, 1),
    "t": Fraction(10**12, 1),
}


@dataclass
class LinearExpr:
    circuit: "Circuit"
    coeffs: dict[int, int] = field(default_factory=dict)
    const: int | None = None

    def __post_init__(self) -> None:
        if self.const is None:
            self.const = self.circuit.spoly.zero
        self.clean()

    def clean(self) -> None:
        zero = self.circuit.spoly.zero
        self.coeffs = {key: value for key, value in self.coeffs.items() if value != zero}
        if self.const is None:
            self.const = zero

    def copy(self) -> "LinearExpr":
        return LinearExpr(self.circuit, dict(self.coeffs), self.const)

    def add(self, other: "LinearExpr") -> "LinearExpr":
        out = dict(self.coeffs)
        for unknown_id, coeff in other.coeffs.items():
            out[unknown_id] = self.circuit.spoly.add(out.get(unknown_id, self.circuit.spoly.zero), coeff)
        return LinearExpr(self.circuit, out, self.circuit.spoly.add(self.const, other.const))

    def neg(self) -> "LinearExpr":
        return LinearExpr(self.circuit, {idx: self.circuit.spoly.neg(coeff) for idx, coeff in self.coeffs.items()}, self.circuit.spoly.neg(self.const))

    def sub(self, other: "LinearExpr") -> "LinearExpr":
        return self.add(other.neg())

    def scale(self, scalar: int) -> "LinearExpr":
        return LinearExpr(
            self.circuit,
            {idx: self.circuit.spoly.mul(scalar, coeff) for idx, coeff in self.coeffs.items()},
            self.circuit.spoly.mul(scalar, self.const),
        )



@dataclass
class Circuit:
    eq_circuit: bool = False
    body: str = "off"
    off_device: set[str] = field(default_factory=set)
    off_device_cap: set[str] = field(default_factory=set)
    zero_cap: set[str] = field(default_factory=set)
    zero_resistance: set[str] = field(default_factory=set)
    zero_inductance: set[str] = field(default_factory=set)
    atoms: AtomTable = field(default_factory=AtomTable)
    expr: ExprArena = field(init=False)
    spoly: SPolyArena = field(init=False)
    symbolic_boundaries: set[str] = field(default_factory=set)
    ground_boundaries: set[str] = field(default_factory=set)
    input_nodes: set[str] = field(default_factory=set)
    device_kinds: dict[str, str] = field(default_factory=dict)
    device_aliases: dict[str, str] = field(default_factory=dict)
    device_currents: dict[str, LinearExpr] = field(default_factory=dict)
    node_set: set[str] = field(default_factory=set)
    node_unknowns: dict[str, int] = field(default_factory=dict)
    branch_unknowns: dict[str, int] = field(default_factory=dict)
    unknowns: list[UnknownInfo] = field(default_factory=list)
    kcl_rows: dict[str, LinearExpr] = field(default_factory=dict)
    constraints: list[LinearExpr] = field(default_factory=list)
    equations: list[LinearExpr] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.expr = ExprArena(self.atoms)
        self.spoly = SPolyArena(self.expr)

    def zero(self) -> LinearExpr:
        return LinearExpr(self)

    def constExpr(self, poly_id: int) -> LinearExpr:
        return LinearExpr(self, {}, poly_id)

    def atomPoly(self, name: str) -> int:
        return self.spoly.fromExpr(self.expr.atom(self.atoms.param(name)))

    def invAtomPoly(self, name: str) -> int:
        return self.spoly.fromExpr(self.expr.atom(self.atoms.invParam(name)))

    def sScaledAtomPoly(self, name: str) -> int:
        return self.spoly.poly((self.expr.zero, self.expr.atom(self.atoms.param(name))))

    def capPoly(self, name: str) -> int:
        if self.isZeroCap(name):
            return self.spoly.zero
        return self.sScaledAtomPoly(name)

    def inductancePoly(self, name: str) -> int:
        if self.isZeroInductance(name):
            return self.spoly.zero
        return self.sScaledAtomPoly(name)

    def dev(self, name: str) -> str:
        key = plain(name)
        return self.device_aliases.get(key, deviceName([name]))

    def isGround(self, node: str) -> bool:
        return nodeName(node) in ground_names

    def isBoundary(self, node: str) -> bool:
        key = nodeName(node)
        return key in self.symbolic_boundaries or key in self.ground_boundaries

    def isKclSkippedNode(self, node: str) -> bool:
        key = nodeName(node)
        return key in ground_names or key in self.symbolic_boundaries or key in self.ground_boundaries

    def nodeVoltage(self, node: str) -> LinearExpr:
        key = nodeName(node)
        self.node_set.add(key)
        if key in ground_names or key in self.ground_boundaries:
            return self.zero()
        if key in self.symbolic_boundaries or key in self.input_nodes:
            return self.constExpr(self.atomPoly(f"v_{{{key}}}"))
        unknown_id = self.node_unknowns.get(key)
        if unknown_id is None:
            unknown_id = len(self.unknowns)
            self.node_unknowns[key] = unknown_id
            self.unknowns.append(UnknownInfo(name=f"v_{{{key}}}", kind="node", node=key))
        return LinearExpr(self, {unknown_id: self.spoly.one}, self.spoly.zero)

    def nodeDiff(self, pos_node: str, neg_node: str) -> LinearExpr:
        return self.nodeVoltage(pos_node).sub(self.nodeVoltage(neg_node))

    def branchCurrentUnknown(self, name: str) -> LinearExpr:
        dev_name = self.dev(name)
        unknown_id = self.branch_unknowns.get(dev_name)
        if unknown_id is None:
            unknown_id = len(self.unknowns)
            self.branch_unknowns[dev_name] = unknown_id
            self.unknowns.append(UnknownInfo(name=f"i_{{{dev_name}}}", kind="branch_current", device=dev_name))
        return LinearExpr(self, {unknown_id: self.spoly.one}, self.spoly.zero)

    def inputCurrent(self, name: str) -> LinearExpr:
        dev_name = self.dev(name)
        return self.constExpr(self.atomPoly(f"i_{{{dev_name}}}"))

    def inputVoltage(self, name: str) -> LinearExpr:
        dev_name = self.dev(name)
        return self.constExpr(self.atomPoly(f"V_{{{dev_name}}}"))

    def addKcl(self, node: str, current: LinearExpr) -> bool:
        key = nodeName(node)
        self.node_set.add(key)
        if self.isKclSkippedNode(node):
            return False
        if not current.coeffs and current.const == self.spoly.zero:
            return False
        old = self.kcl_rows.get(key, self.zero())
        self.kcl_rows[key] = old.add(current)
        return True

    def addConstraint(self, equation: LinearExpr) -> None:
        if equation.coeffs or equation.const != self.spoly.zero:
            self.constraints.append(equation)

    def finish(self) -> "Circuit":
        self.equations = [self.kcl_rows[key] for key in sorted(self.kcl_rows)] + self.constraints
        return self

    def recordDeviceCurrent(self, name: str, current: LinearExpr) -> None:
        self.device_currents[self.dev(name)] = current

    def isConductanceOff(self, name: str) -> bool:
        return self.dev(name) in self.off_device

    def isCapacitanceOff(self, name: str) -> bool:
        return self.dev(name) in self.off_device_cap

    def isZeroCap(self, name: str) -> bool:
        return plain(name) in self.zero_cap

    def isZeroResistance(self, name: str) -> bool:
        key = plain(name)
        return key in self.zero_resistance or f"inv_{key}" in self.zero_resistance or f"1_{key}" in self.zero_resistance

    def isZeroInductance(self, name: str) -> bool:
        return plain(name) in self.zero_inductance


def buildCircuit(
    records: list[list[str]],
    eq_circuit: bool = False,
    body: str = "off",
    off_device=None,
    off_device_cap=None,
    zero_cap=None,
    zero_resistance=None,
    zero_inductance=None,
) -> Circuit:
    if not isinstance(eq_circuit, bool):
        raise ValueError("eq_circuit must be True or False")
    body = str(body).lower()
    if body not in {"off", "on"}:
        raise ValueError("body must be 'off' or 'on'")
    errors = [" ".join(tokens[1:]) for tokens in records if tokens and tokens[0] == "__error__"]
    if errors:
        raise ValueError("invalid netlist: " + "; ".join(errors))
    symbolic_boundaries = {nodeName(node) for tokens in records if tokens and tokens[0] == "__boundary_symbolic__" for node in tokens[1:]}
    ground_boundaries = {nodeName(node) for tokens in records if tokens and tokens[0] == "__boundary_ground__" for node in tokens[1:]}
    boundary_markers = {"__boundary_symbolic__", "__boundary_ground__"}
    primitive_records = [tokens for tokens in records if tokens and tokens[0] not in boundary_markers]
    zero_resistance_items = optionItems(zero_resistance, "zero_resistance")
    circuit = Circuit(
        eq_circuit=eq_circuit,
        body=body,
        symbolic_boundaries=symbolic_boundaries,
        ground_boundaries=ground_boundaries,
        zero_cap=variableOptionSet(zero_cap, "zero_cap"),
        zero_resistance=variableOptionSet(zero_resistance_items, "zero_resistance"),
        zero_inductance=variableOptionSet(zero_inductance, "zero_inductance"),
    )
    bindDevices(circuit, primitive_records)
    compact_devices = {circuit.dev(tokens[0]) for tokens in primitive_records if instKind(tokens) in {"m", "q", "d"}}
    circuit.off_device = deviceOptionSet(circuit, off_device, compact_devices, "off_device")
    circuit.off_device_cap = deviceOptionSet(circuit, off_device_cap, compact_devices, "off_device_cap")
    markInputVoltageNodes(circuit, primitive_records)
    for tokens in primitive_records:
        stamper = stamp_map.get(instKind(tokens))
        if stamper is not None:
            stamper(circuit, tokens)
    return circuit.finish()


def bindDevices(circuit: Circuit, records: list[list[str]]) -> None:
    for tokens in records:
        canonical = deviceName(tokens)
        kind = instKind(tokens)
        circuit.device_kinds[canonical] = kind
        for alias in deviceAliases(tokens[0], canonical):
            circuit.device_aliases[plain(alias)] = canonical
        for idx in nodeIndices(kind):
            if idx < len(tokens):
                circuit.node_set.add(nodeName(tokens[idx]))


def optionItems(value, option_name: str) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    try:
        items = list(value)
    except TypeError as exc:
        raise ValueError(f"{option_name} must be None, a string, or an iterable of strings") from exc
    if not all(isinstance(item, str) for item in items):
        raise ValueError(f"{option_name} entries must be strings")
    return items


def deviceOptionSet(circuit: Circuit, value, compact_devices: set[str], option_name: str) -> set[str]:
    devices = {circuit.dev(item) for item in optionItems(value, option_name)}
    unknown = sorted(device for device in devices if device not in compact_devices)
    if unknown:
        raise ValueError(f"{option_name}: unknown compact device(s): " + ", ".join(unknown))
    return devices


def variableOptionSet(value, option_name: str) -> set[str]:
    return {plain(item) for item in optionItems(value, option_name) if plain(item)}


def markInputVoltageNodes(circuit: Circuit, records: list[list[str]]) -> None:
    for tokens in records:
        if instKind(tokens) != "v" or isZeroAc(tokens[3:]):
            continue
        pos_node, neg_node = tokens[1], tokens[2]
        if circuit.isGround(neg_node) and not circuit.isBoundary(pos_node):
            circuit.input_nodes.add(nodeName(pos_node))
        elif circuit.isGround(pos_node) and not circuit.isBoundary(neg_node):
            circuit.input_nodes.add(nodeName(neg_node))


def stampVoltageSource(circuit: Circuit, tokens: list[str]) -> None:
    name, pos_node, neg_node = tokens[:3]
    current = circuit.branchCurrentUnknown(name)
    if stampCurrent(circuit, pos_node, neg_node, current):
        circuit.recordDeviceCurrent(name, current)
    if isZeroAc(tokens[3:]):
        circuit.addConstraint(circuit.nodeDiff(pos_node, neg_node))
        return
    if groundReferencedInputNode(circuit, pos_node, neg_node) is not None:
        return
    circuit.addConstraint(circuit.nodeDiff(pos_node, neg_node).sub(circuit.inputVoltage(name)))


def groundReferencedInputNode(circuit: Circuit, pos_node: str, neg_node: str) -> str | None:
    if circuit.isGround(neg_node) and not circuit.isBoundary(pos_node):
        return pos_node
    if circuit.isGround(pos_node) and not circuit.isBoundary(neg_node):
        return neg_node
    return None


def stampCurrentSource(circuit: Circuit, tokens: list[str]) -> None:
    name, pos_node, neg_node = tokens[:3]
    if isZeroAc(tokens[3:]):
        circuit.recordDeviceCurrent(name, circuit.zero())
        return
    current = circuit.inputCurrent(name)
    circuit.recordDeviceCurrent(name, current)
    stampCurrent(circuit, pos_node, neg_node, current)


def stampVcvSource(circuit: Circuit, tokens: list[str]) -> None:
    name, pos_node, neg_node, ctrl_pos, ctrl_neg, gain_token = tokens[:6]
    current = circuit.branchCurrentUnknown(name)
    stampCurrent(circuit, pos_node, neg_node, current)
    gain = parseValue(circuit, gain_token)
    circuit.addConstraint(circuit.nodeDiff(pos_node, neg_node).sub(circuit.nodeDiff(ctrl_pos, ctrl_neg).scale(gain)))


def stampVccSource(circuit: Circuit, tokens: list[str]) -> None:
    name, pos_node, neg_node, ctrl_pos, ctrl_neg, gain_token = tokens[:6]
    gain = parseValue(circuit, gain_token)
    current = circuit.nodeDiff(ctrl_pos, ctrl_neg).scale(gain)
    stampCurrent(circuit, pos_node, neg_node, current)
    circuit.recordDeviceCurrent(name, current)


def stampCccSource(circuit: Circuit, tokens: list[str]) -> None:
    name, pos_node, neg_node, ctrl_source, gain_token = tokens[:5]
    gain = parseValue(circuit, gain_token)
    ctrl_current = circuit.branchCurrentUnknown(ctrl_source)
    current = ctrl_current.scale(gain)
    stampCurrent(circuit, pos_node, neg_node, current)
    circuit.recordDeviceCurrent(name, current)


def stampCcvSource(circuit: Circuit, tokens: list[str]) -> None:
    name, pos_node, neg_node, ctrl_source, gain_token = tokens[:5]
    current = circuit.branchCurrentUnknown(name)
    ctrl_current = circuit.branchCurrentUnknown(ctrl_source)
    stampCurrent(circuit, pos_node, neg_node, current)
    gain = parseValue(circuit, gain_token)
    circuit.addConstraint(circuit.nodeDiff(pos_node, neg_node).sub(ctrl_current.scale(gain)))
    circuit.recordDeviceCurrent(name, current)


def stampResistor(circuit: Circuit, tokens: list[str]) -> None:
    name, pos_node, neg_node = tokens[:3]
    stampResistance(circuit, name, pos_node, neg_node, passive("R", name))


def stampCapacitor(circuit: Circuit, tokens: list[str]) -> None:
    name, pos_node, neg_node = tokens[:3]
    stampAdmittance(circuit, name, pos_node, neg_node, circuit.capPoly(passive("C", name)))


def stampInductor(circuit: Circuit, tokens: list[str]) -> None:
    name, pos_node, neg_node = tokens[:3]
    current = circuit.branchCurrentUnknown(name)
    if stampCurrent(circuit, pos_node, neg_node, current):
        circuit.recordDeviceCurrent(name, current)
    voltage = circuit.nodeDiff(pos_node, neg_node)
    s_l = circuit.inductancePoly(passive("L", name))
    circuit.addConstraint(voltage.sub(current.scale(s_l)))


def stampDiode(circuit: Circuit, tokens: list[str]) -> None:
    name, pos_node, neg_node = tokens[:3]
    inst = circuit.dev(name)
    current = circuit.zero()
    voltage = circuit.nodeDiff(pos_node, neg_node)
    if not circuit.isConductanceOff(inst):
        current = current.add(voltage.scale(circuit.atomPoly(f"g_{{{inst}}}")))
    if not circuit.isCapacitanceOff(inst):
        current = current.add(voltage.scale(circuit.capPoly(f"C_{{{inst}}}")))
    circuit.recordDeviceCurrent(name, current)
    stampCurrent(circuit, pos_node, neg_node, current)


def stampMos(circuit: Circuit, tokens: list[str]) -> None:
    name, drain, gate, source, body = tokens[:5]
    nodes = {"g": gate, "d": drain, "s": source, "b": body}
    if circuit.eq_circuit:
        stampMosEq(circuit, name, nodes)
    else:
        stampMosBsim(circuit, name, nodes)


def stampBjt(circuit: Circuit, tokens: list[str]) -> None:
    name, collector, base, emitter = tokens[:4]
    nodes = {"c": collector, "b": base, "e": emitter}
    if circuit.eq_circuit:
        stampBjtEq(circuit, name, nodes)
    else:
        stampBjtDifferential(circuit, name, nodes)


def stampCurrent(circuit: Circuit, pos_node: str, neg_node: str, current: LinearExpr) -> bool:
    pos_used = circuit.addKcl(pos_node, current)
    neg_used = circuit.addKcl(neg_node, current.neg())
    return pos_used or neg_used


def stampAdmittance(circuit: Circuit, name: str, pos_node: str, neg_node: str, admittance: int) -> bool:
    if admittance == circuit.spoly.zero:
        circuit.recordDeviceCurrent(name, circuit.zero())
        return False
    voltage = circuit.nodeDiff(pos_node, neg_node)
    current = voltage.scale(admittance)
    circuit.recordDeviceCurrent(name, current)
    return stampCurrent(circuit, pos_node, neg_node, current)


def stampShort(circuit: Circuit, name: str, pos_node: str, neg_node: str) -> None:
    current = circuit.branchCurrentUnknown(name)
    if stampCurrent(circuit, pos_node, neg_node, current):
        circuit.recordDeviceCurrent(name, current)
    circuit.addConstraint(circuit.nodeDiff(pos_node, neg_node))


def stampResistance(circuit: Circuit, name: str, pos_node: str, neg_node: str, resistance_name: str) -> None:
    if circuit.isZeroResistance(resistance_name):
        stampShort(circuit, name, pos_node, neg_node)
        return
    stampAdmittance(circuit, name, pos_node, neg_node, circuit.invAtomPoly(resistance_name))


def stampControlledCurrent(circuit: Circuit, pos_node: str, neg_node: str, coeff_name: str, ctrl_pos: str, ctrl_neg: str) -> bool:
    coeff = circuit.atomPoly(coeff_name)
    if coeff == circuit.spoly.zero:
        return False
    current = circuit.nodeDiff(ctrl_pos, ctrl_neg).scale(coeff)
    return stampCurrent(circuit, pos_node, neg_node, current)


def stampMosBsim(circuit: Circuit, name: str, nodes: dict[str, str]) -> None:
    inst = circuit.dev(name)
    rows = mos_rows if circuit.body == "on" else ("g", "d", "s")
    for row in rows:
        if circuit.isKclSkippedNode(nodes[row]):
            continue
        total = circuit.zero()
        if not circuit.isConductanceOff(inst):
            coeff = circuit.atomPoly(f"g_{{{inst},gs,{row}}}")
            if coeff != circuit.spoly.zero:
                total = total.add(circuit.nodeDiff(nodes["g"], nodes["s"]).scale(coeff))
            coeff = circuit.atomPoly(f"g_{{{inst},ds,{row}}}")
            if coeff != circuit.spoly.zero:
                total = total.add(circuit.nodeDiff(nodes["d"], nodes["s"]).scale(coeff))
            if circuit.body == "on":
                coeff = circuit.atomPoly(f"g_{{{inst},bs,{row}}}")
                if coeff != circuit.spoly.zero:
                    total = total.add(circuit.nodeDiff(nodes["b"], nodes["s"]).scale(coeff))
        if not circuit.isCapacitanceOff(inst):
            for col in rows:
                coeff = circuit.capPoly(f"C_{{{inst},{col},{row}}}")
                if coeff != circuit.spoly.zero:
                    total = total.add(circuit.nodeVoltage(nodes[col]).scale(coeff))
        circuit.addKcl(nodes[row], total)


def stampBjtDifferential(circuit: Circuit, name: str, nodes: dict[str, str]) -> None:
    inst = circuit.dev(name)
    for row in bjt_rows:
        if circuit.isKclSkippedNode(nodes[row]):
            continue
        total = circuit.zero()
        if not circuit.isConductanceOff(inst):
            coeff = circuit.atomPoly(f"g_{{{inst},be,{row}}}")
            if coeff != circuit.spoly.zero:
                total = total.add(circuit.nodeDiff(nodes["b"], nodes["e"]).scale(coeff))
            coeff = circuit.atomPoly(f"g_{{{inst},ce,{row}}}")
            if coeff != circuit.spoly.zero:
                total = total.add(circuit.nodeDiff(nodes["c"], nodes["e"]).scale(coeff))
        if not circuit.isCapacitanceOff(inst):
            for col in bjt_rows:
                coeff = circuit.capPoly(f"C_{{{inst},{col},{row}}}")
                if coeff != circuit.spoly.zero:
                    total = total.add(circuit.nodeVoltage(nodes[col]).scale(coeff))
        circuit.addKcl(nodes[row], total)


def stampMosEq(circuit: Circuit, name: str, nodes: dict[str, str]) -> None:
    inst = circuit.dev(name)
    if not circuit.isConductanceOff(inst):
        stampControlledCurrent(circuit, nodes["d"], nodes["s"], f"g_{{{inst}}}", nodes["g"], nodes["s"])
        if circuit.body == "on":
            stampControlledCurrent(circuit, nodes["d"], nodes["s"], f"g_{{{inst},mb}}", nodes["b"], nodes["s"])
        stampResistance(circuit, f"{inst}_r", nodes["d"], nodes["s"], f"r_{{{inst}}}")
    stampCapPairs(circuit, inst, nodes, mos_body_cap_pairs if circuit.body == "on" else mos_cap_pairs)


def stampBjtEq(circuit: Circuit, name: str, nodes: dict[str, str]) -> None:
    inst = circuit.dev(name)
    if not circuit.isConductanceOff(inst):
        stampControlledCurrent(circuit, nodes["c"], nodes["e"], f"g_{{{inst}}}", nodes["b"], nodes["e"])
        stampResistance(circuit, f"{inst}_r_o", nodes["c"], nodes["e"], f"r_{{{inst},o}}")
        stampResistance(circuit, f"{inst}_pi", nodes["b"], nodes["e"], f"r_{{{inst},pi}}")
    stampCapPairs(circuit, inst, nodes, bjt_cap_pairs)


def stampCapPairs(circuit: Circuit, inst: str, nodes: dict[str, str], pairs: dict[str, tuple[str, str]]) -> None:
    if circuit.isCapacitanceOff(inst):
        return
    for label, (pos_key, neg_key) in pairs.items():
        pos_node, neg_node = nodes[pos_key], nodes[neg_key]
        if circuit.isKclSkippedNode(pos_node) and circuit.isKclSkippedNode(neg_node):
            continue
        stampAdmittance(circuit, f"{inst}_{label}", pos_node, neg_node, circuit.capPoly(f"C_{{{inst},{label}}}"))


def parseValue(circuit: Circuit, token: str) -> int:
    text = token.split("=", 1)[-1].strip().strip("{}")
    frac = parseNumber(text)
    if frac is not None:
        return circuit.spoly.fromExpr(circuit.expr.rational(frac))
    return circuit.atomPoly(text)


def parseNumber(text: str) -> Fraction | None:
    raw = text.strip().lower()
    if not raw:
        return Fraction(0)
    suffix = ""
    number = raw
    for candidate in sorted(unit_scale, key=len, reverse=True):
        if raw.endswith(candidate) and len(raw) > len(candidate):
            suffix = candidate
            number = raw[: -len(candidate)]
            break
    try:
        value = Fraction(number)
    except ValueError:
        try:
            value = Fraction(float(number)).limit_denominator(10**12)
        except ValueError:
            return None
    return value * unit_scale.get(suffix, Fraction(1))


def acValue(tokens: list[str]) -> str:
    lowered = [token.lower() for token in tokens]
    for idx, token in enumerate(lowered):
        if token == "ac" and idx + 1 < len(tokens):
            return tokens[idx + 1]
        if token.startswith("ac="):
            return tokens[idx][3:].split(",", 1)[0]
    return "0"


def isZeroAc(tokens: list[str]) -> bool:
    text = acValue(tokens).strip().strip("{}")
    value = parseNumber(text)
    if value is not None:
        return value == 0
    return text in {"", "0"}


def passive(prefix: str, name: str) -> str:
    return f"{prefix}_{{{suffixName(rawSuffix(name))}}}"


def rawSuffix(name: str) -> str:
    clean = plain(name)
    return clean[1:] if len(clean) > 1 else clean


def suffixName(text: str) -> str:
    return text if len(text) == 1 else text.lower()


def instKind(tokens: list[str]) -> str:
    low_name = tokens[0].lower()
    if low_name.startswith(("mn", "mp", "nm", "pm", "m")):
        return "m"
    if low_name.startswith(("qn", "qp", "nq", "pq", "q")):
        return "q"
    return low_name[:1]


def deviceName(tokens: list[str]) -> str:
    name = plain(tokens[0])
    low_name = name.lower()
    if low_name.startswith(("m", "q", "d")):
        return low_name
    return familyName(name)


def familyName(name: str) -> str:
    return name.upper() if len(name) <= 1 else name[0].upper() + "_" + suffixName(name[1:])


def mosSuffix(name: str, low_name: str) -> str:
    return name[2:].lower() if low_name.startswith(("mn", "mp", "nm", "pm")) else name[1:].lower()


def bjtSuffix(name: str, low_name: str) -> str:
    return name[2:].lower() if low_name.startswith(("qn", "qp", "nq", "pq")) else name[1:].lower()


def deviceAliases(raw: str, canonical: str) -> set[str]:
    aliases = {canonical, plain(raw)}
    low_name = plain(raw).lower()
    if low_name.startswith(("mn", "nm")):
        aliases.add("M_n" + mosSuffix(plain(raw), low_name))
    elif low_name.startswith(("mp", "pm")):
        aliases.add("M_p" + mosSuffix(plain(raw), low_name))
    elif low_name.startswith(("qn", "nq")):
        aliases.add("Q_n" + bjtSuffix(plain(raw), low_name))
    elif low_name.startswith(("qp", "pq")):
        aliases.add("Q_p" + bjtSuffix(plain(raw), low_name))
    return aliases


def nodeIndices(kind: str) -> tuple[int, ...]:
    return {
        "v": (1, 2),
        "i": (1, 2),
        "r": (1, 2),
        "c": (1, 2),
        "l": (1, 2),
        "d": (1, 2),
        "e": (1, 2, 3, 4),
        "g": (1, 2, 3, 4),
        "f": (1, 2),
        "h": (1, 2),
        "m": (1, 2, 3, 4),
        "q": (1, 2, 3),
    }.get(kind, ())


stamp_map = {
    "v": stampVoltageSource,
    "i": stampCurrentSource,
    "e": stampVcvSource,
    "g": stampVccSource,
    "f": stampCccSource,
    "h": stampCcvSource,
    "r": stampResistor,
    "c": stampCapacitor,
    "l": stampInductor,
    "d": stampDiode,
    "m": stampMos,
    "q": stampBjt,
}
