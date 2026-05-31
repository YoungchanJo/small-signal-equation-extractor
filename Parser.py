from pathlib import Path
import re


INSTANCE_PREFIXES = ("mn", "mp", "nm", "pm", "qn", "qp", "nq", "pq", "m", "q", "v", "i", "e", "g", "f", "h", "r", "c", "l", "d")
MIN_FIELDS = {"v": 4, "i": 4, "e": 6, "g": 6, "f": 5, "h": 5, "r": 4, "c": 4, "l": 4, "d": 4, "m": 6, "q": 5}
NODE_IDX = {
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
}
SKIP_BLOCK_END = {".control": ".endc", ".data": ".enddata"}
SUBCKT_BOUNDARY_MARKER = {"ground": "__boundary_ground__", "symbolic": "__boundary_symbolic__"}


def parseLines(
    text: str,
    subckt: str | None = None,
    subckt_boundary: str = "ground",
) -> list[list[str]]:
    if subckt is not None and not isinstance(subckt, str):
        return errorRecord("subckt must be a string or None")
    if subckt_boundary not in SUBCKT_BOUNDARY_MARKER:
        return errorRecord("subckt_boundary must be 'ground' or 'symbolic'")
    main_records, sub_defs, x_calls = scanLines(text)
    if subckt is None:
        return main_records

    key = subckt.lower()
    if key not in sub_defs:
        return errorRecord(f"subckt not found: {subckt}")
    if len(sub_defs[key]) != 1:
        return errorRecord(f"unsupported duplicate subckt: {subckt}")

    pins, sub_records = sub_defs[key][0]
    mapping = selectedSubcktPinMap(pins, x_calls, key)
    if not mapping:
        return errorRecord(f"subckt instance not found: {subckt}")

    mapped = [mapNodes(tokens, mapping) for tokens in sub_records]
    duplicates = duplicateDevices(main_records, mapped)
    if duplicates:
        return errorRecord("duplicate device names between main and subckt: " + ", ".join(duplicates))

    boundary_nodes = skippedSubcktBoundaryNodes(main_records, x_calls, key, set(mapping.values()))
    records = ([[SUBCKT_BOUNDARY_MARKER[subckt_boundary], *boundary_nodes]] if boundary_nodes else []) + main_records + mapped
    return records


def parseFile(
    path: str,
    subckt: str | None = None,
    subckt_boundary: str = "ground",
) -> list[list[str]]:
    file_path = Path(path)
    if not file_path.is_file():
        return errorRecord(f"file not found: {path}")
    return parseLines(file_path.read_text(), subckt=subckt, subckt_boundary=subckt_boundary)


def scanLines(
    text: str,
) -> tuple[list[list[str]], dict[str, list[tuple[list[str], list[list[str]]]]], list[list[str]]]:
    main_records: list[list[str]] = []
    sub_defs: dict[str, list[tuple[list[str], list[list[str]]]]] = {}
    x_calls: list[list[str]] = []
    sub_name, sub_pins, sub_records, skip_until = "", [], [], ""
    for raw_line in text.splitlines():
        line = lineCode(raw_line).strip()
        if not line:
            continue
        tokens = line.split()
        head = tokens[0].lower()

        if sub_name:
            if isSubcktEnd(head):
                sub_defs.setdefault(sub_name.lower(), []).append((sub_pins, sub_records))
                sub_name, sub_pins, sub_records = "", [], []
            else:
                sub_records.append(parseInstanceRecord(tokens))
            continue

        if skip_until:
            skip_until = "" if head == skip_until else skip_until
            continue
        if head == ".end":
            break
        if head == ".subckt" and len(tokens) >= 2:
            sub_name, sub_pins, sub_records = tokens[1], tokens[2:], []
            continue
        if head in SKIP_BLOCK_END:
            skip_until = SKIP_BLOCK_END[head]
            continue
        if head.startswith("."):
            continue
        if head.startswith("x"):
            x_calls.append(tokens)
            continue
        main_records.append(parseInstanceRecord(tokens))

    return main_records, sub_defs, x_calls


def parseInstanceRecord(tokens: list[str]) -> list[str]:
    error = instanceError(tokens)
    if error is None and instanceKind(tokens[0]):
        return tokens
    if error is not None:
        return ["__error__", error]
    return ["__error__", f"unsupported line: {tokens[0]}"]


def instanceError(tokens: list[str]) -> str | None:
    kind = instanceKind(tokens[0])
    if not kind:
        return None
    if len(tokens) < MIN_FIELDS[kind] or (kind in {"r", "c", "l"} and len(tokens) != MIN_FIELDS[kind]):
        return f"{tokens[0]}: invalid {kind.upper()} instance"
    for idx in NODE_IDX.get(kind, ()):
        if idx >= len(tokens) or not validName(tokens[idx]):
            return f"{tokens[0]}: invalid node name '{tokens[idx] if idx < len(tokens) else ''}'"
    if not validName(tokens[0]):
        return f"{tokens[0]}: invalid device name"
    return None


def validName(text: str) -> bool:
    return bool(cleanName(text))


def isSubcktEnd(head: str) -> bool:
    return head == ".end" or head.startswith(".ends")


def selectedSubcktPinMap(pins: list[str], x_calls: list[list[str]], subckt_key: str) -> dict[str, str]:
    for tokens in x_calls:
        if tokens[-1].lower() == subckt_key and len(tokens) >= len(pins) + 2:
            return dict(zip((pin.lower() for pin in pins), tokens[1 : 1 + len(pins)]))
    return {}


def skippedSubcktBoundaryNodes(
    main_records: list[list[str]],
    x_calls: list[list[str]],
    selected_key: str,
    selected_nodes: set[str],
) -> list[str]:
    main_nodes = {
        tokens[idx].lower()
        for tokens in main_records
        if tokens and tokens[0] != "__error__"
        for idx in NODE_IDX.get(instanceKind(tokens[0]), ())
        if idx < len(tokens)
    }
    selected_lower = {node.lower() for node in selected_nodes}
    boundary = set()
    for tokens in x_calls:
        if not tokens or tokens[-1].lower() == selected_key:
            continue
        boundary |= {
            node
            for node in tokens[1:-1]
            if node.lower() in selected_lower and node.lower() not in main_nodes
        }
    return sorted(boundary, key=str.lower)


def mapNodes(tokens: list[str], mapping: dict[str, str]) -> list[str]:
    if not tokens or tokens[0] == "__error__":
        return tokens.copy()
    out = tokens.copy()
    for idx in NODE_IDX.get(instanceKind(tokens[0]), ()):
        if idx < len(out) and out[idx].lower() in mapping:
            out[idx] = mapping[out[idx].lower()]
    return out


def duplicateDevices(main_records: list[list[str]], sub_records: list[list[str]]) -> list[str]:
    main = {deviceKey(tokens) for tokens in main_records if tokens and tokens[0] != "__error__"}
    sub = {deviceKey(tokens) for tokens in sub_records if tokens and tokens[0] != "__error__"}
    return sorted(main & sub)


def deviceKey(tokens: list[str]) -> str:
    return cleanName(tokens[0]).lower()


def cleanName(text: str) -> str:
    return re.sub(r"[^0-9A-Za-z]+", "_", str(text)).strip("_")


def errorRecord(message: str) -> list[list[str]]:
    return [["__error__", message]]


def instanceKind(name: str) -> str:
    name = name.lower()
    if name.startswith(("mn", "mp", "nm", "pm", "m")):
        return "m"
    if name.startswith(("qn", "qp", "nq", "pq", "q")):
        return "q"
    return name[:1] if name.startswith(INSTANCE_PREFIXES) else ""


def lineCode(line: str) -> str:
    if not line.strip() or line.lstrip().startswith("*"):
        return ""
    cuts = [len(line)] + [idx for mark in ("//", ";", "$") if (idx := line.find(mark)) >= 0]
    return line[: min(cuts)]
