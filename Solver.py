from __future__ import annotations

import re

from SolverKernel import KernelRequestError, Response, SolverKernel, splitTopLevelSlash, nodeName


probe_pattern = re.compile(r"^\s*([VI])\s*\(\s*([^()]+?)\s*\)\s*$", re.IGNORECASE)


class Solver:
    def __init__(self, circuit):
        self.circuit = circuit
        self.kernel = SolverKernel(circuit)

    def response(self, request: str) -> Response:
        if not isinstance(request, str):
            raise TypeError("request must be a string")
        parts = splitTopLevelSlash(request.strip())
        if parts is None:
            numerator_text = request.strip()
            denominator_probe = self.circuit.constExpr(self.circuit.spoly.one)
        else:
            numerator_text, denominator_text = parts
            denominator_probe = self.parseProbe(denominator_text)
        numerator_probe = self.parseProbe(numerator_text)
        return self.kernel.responseFromProbes(request.strip(), numerator_probe, denominator_probe)

    def parseProbe(self, text: str):
        match = probe_pattern.match(text)
        if not match:
            raise KernelRequestError(f"invalid probe syntax: {text}")
        kind = match.group(1).upper()
        name = match.group(2).strip()
        if kind == "V":
            return self.nodeVoltageRequest(name)
        return self.deviceCurrentRequest(name)

    def nodeVoltageRequest(self, node: str):
        key = nodeName(node)
        if key not in self.circuit.node_set and key not in {"0", "gnd"}:
            raise KernelRequestError(f"V({node}) refers to an unknown netlist node")
        return self.circuit.nodeVoltage(node)

    def deviceCurrentRequest(self, name: str):
        dev_name = self.circuit.dev(name)
        kind = self.circuit.device_kinds.get(dev_name)
        if kind is None:
            raise KernelRequestError(f"I({name}) refers to an unknown netlist device")
        if kind not in {"v", "i", "r", "c", "l", "d"}:
            raise KernelRequestError(f"I({name}) is defined only for devices with a public two-terminal current convention")
        current = self.circuit.device_currents.get(dev_name)
        if current is None:
            raise KernelRequestError(f"I({name}) has no registered two-terminal current")
        return current
