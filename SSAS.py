from __future__ import annotations

import inspect
from os import PathLike
from pathlib import Path
import re
from typing import Literal

from Nodal import buildCircuit
from Parser import parseFile, parseLines
from Solver import Solver
from SolverKernel import Response, plain


DISPLAY_MATH_BLOCK = re.compile(r"\\\[(.*?)\\\]", re.DOTALL)
SourceMode = Literal["file", "inline"]


class _DisplayMath:
    def __init__(self, body: str):
        self.body = body.strip()

    def _repr_latex_(self) -> str:
        return "$\\displaystyle " + self.body + "$"

    def _repr_markdown_(self) -> str:
        return "\\[\n" + self.body + "\n\\]"

    def __repr__(self) -> str:
        return self._repr_markdown_()


class SSAS:
    """User-facing wrapper for symbolic small-signal extraction.

    The public workflow is solve(...). Each call parses the selected netlist
    source, solves one request, writes one Markdown file, and optionally displays
    the Markdown in an IPython notebook. solve(...) intentionally returns None.
    """

    def __init__(
        self,
        output_dir: str | PathLike[str] = "tex_outputs",
        *,
        absolute_route: bool = False,
    ):
        self.absolute_route = absolute_route
        self.route_dir = self._caller_dir()
        self.output_dir = self._route_path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.circuit = None

    def solve(
        self,
        netlist: str | PathLike[str],
        request: str,
        output_name: str | PathLike[str] | None = None,
        *,
        source: SourceMode,
        eq_circuit: bool = False,
        body: str = "off",
        subckt: str | None = None,
        subckt_boundary: str = "ground",
        off_device=None,
        off_device_cap=None,
        zero_cap=None,
        zero_resistance=None,
        zero_inductance=None,
        equal_passive=None,
        zero_pole: bool = False,
        fraction_symbol_threshold: int = 10,
        show: bool = False,
    ) -> None:
        """Solve one request and save one Markdown file.

        source must be explicit:
        - source="inline": netlist is SPICE-like netlist text.
        - source="file": netlist is a path resolved by the SSAS route rule.
        """
        solver = self._make_solver(
            netlist,
            source=source,
            eq_circuit=eq_circuit,
            body=body,
            subckt=subckt,
            subckt_boundary=subckt_boundary,
            off_device=off_device,
            off_device_cap=off_device_cap,
            zero_cap=zero_cap,
            zero_resistance=zero_resistance,
            zero_inductance=zero_inductance,
        )
        response = self._with_equal_passive(solver.response(request), equal_passive)
        markdown = response.kernel.renderTex(
            response,
            zero_pole=zero_pole,
            fraction_symbol_threshold=fraction_symbol_threshold,
        )
        output_path = self._output_path(output_name, request)
        output_path.write_text(markdown)
        if show:
            self._display_markdown(markdown)

    def _make_solver(
        self,
        netlist: str | PathLike[str],
        *,
        source: SourceMode,
        eq_circuit: bool,
        body: str,
        subckt: str | None,
        subckt_boundary: str,
        off_device,
        off_device_cap,
        zero_cap,
        zero_resistance,
        zero_inductance,
    ) -> Solver:
        if source == "file":
            records = parseFile(str(self._route_path(netlist)), subckt=subckt, subckt_boundary=subckt_boundary)
        elif source == "inline":
            records = parseLines(str(netlist), subckt=subckt, subckt_boundary=subckt_boundary)
        else:
            raise ValueError("source must be 'file' or 'inline'")
        self.circuit = buildCircuit(
            records,
            eq_circuit=eq_circuit,
            body=body,
            off_device=off_device,
            off_device_cap=off_device_cap,
            zero_cap=zero_cap,
            zero_resistance=zero_resistance,
            zero_inductance=zero_inductance,
        )
        return Solver(self.circuit)

    def _with_equal_passive(self, response: Response, equal_passive) -> Response:
        response.passive_value_aliases = self._passive_value_aliases(equal_passive)
        return response

    def _passive_value_aliases(self, value, option_name: str = "equal_passive") -> dict[str, str]:
        if self.circuit is None:
            raise ValueError("no circuit is available")
        aliases: dict[str, str] = {}
        for alias, devices in self._equal_passive_groups(value, option_name):
            if not isinstance(alias, str) or not plain(alias):
                raise ValueError(f"{option_name}: alias must be a non-empty string")
            device_items = self._option_items(devices, f"{option_name} device list")
            if not device_items:
                raise ValueError(f"{option_name}: equal_device list must not be empty")

            resolved: list[tuple[str, str]] = []
            invalid: list[str] = []
            for item in device_items:
                device = self.circuit.dev(item)
                kind = self.circuit.device_kinds.get(device)
                if kind not in {"r", "c", "l"}:
                    invalid.append(item)
                    continue
                resolved.append((device, kind))
            if invalid:
                raise ValueError(f"{option_name}: unknown or non-passive device(s): " + ", ".join(invalid))

            kinds = {kind for _, kind in resolved}
            if len(kinds) != 1:
                raise ValueError(f"{option_name}: one alias group cannot mix R, C, and L devices")
            target = self._passive_alias_symbol(alias, next(iter(kinds)), option_name)

            for device, _ in resolved:
                source = plain(self._passive_symbol_from_canonical(device))
                existing = aliases.get(source)
                if existing is not None and existing != target:
                    raise ValueError(f"{option_name}: passive device {device} is assigned to multiple aliases")
                aliases[source] = target
        return aliases

    def _passive_alias_symbol(self, alias: str, kind: str, option_name: str) -> str:
        assert self.circuit is not None
        device = self.circuit.dev(alias)
        if device in self.circuit.device_kinds:
            alias_kind = self.circuit.device_kinds[device]
            if alias_kind not in {"r", "c", "l"}:
                raise ValueError(f"{option_name}: alias cannot be a non-passive device")
            if alias_kind != kind:
                raise ValueError(f"{option_name}: alias passive type must match equal_device type")
            return self._passive_symbol_from_canonical(device)
        return alias

    @staticmethod
    def _equal_passive_groups(value, option_name: str) -> list[tuple[object, object]]:
        if value is None:
            return []
        if isinstance(value, dict):
            if "alias" in value and "equal_device" in value:
                return [(value["alias"], value["equal_device"])]
            raise ValueError(f"{option_name} dict groups must contain alias and equal_device")
        if isinstance(value, (list, tuple)) and len(value) == 2 and isinstance(value[0], str):
            return [(value[0], value[1])]
        try:
            groups = list(value)
        except TypeError as exc:
            raise ValueError(f"{option_name} must be None, an (alias, equal_device) pair, or an iterable of pairs") from exc

        parsed: list[tuple[object, object]] = []
        for group in groups:
            if isinstance(group, dict):
                if "alias" not in group or "equal_device" not in group:
                    raise ValueError(f"{option_name} dict groups must contain alias and equal_device")
                parsed.append((group["alias"], group["equal_device"]))
            elif isinstance(group, (list, tuple)) and len(group) == 2:
                parsed.append((group[0], group[1]))
            else:
                raise ValueError(f"{option_name} groups must be (alias, equal_device) pairs")
        return parsed

    @staticmethod
    def _option_items(value, option_name: str) -> list[str]:
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

    @staticmethod
    def _passive_symbol_from_canonical(canonical: str) -> str:
        key = plain(canonical)
        prefix = key[:1].upper()
        if "_" in key:
            suffix = key.split("_", 1)[1]
        else:
            suffix = key[1:] or key
        return f"{prefix}_{{{suffix}}}"

    def _output_path(
        self,
        output_name: str | PathLike[str] | None,
        request: str,
    ) -> Path:
        if output_name is None:
            output_name = self._safe_stem(request) + ".md"
        path = Path(output_name)
        if not path.suffix:
            path = path.with_suffix(".md")
        if not path.is_absolute():
            path = self.output_dir / path
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    @staticmethod
    def _display_markdown(text: str) -> None:
        Markdown = SSAS._markdown_type()
        if Markdown is None:
            print(text)
            return
        cursor = 0
        rendered_math = False
        for match in DISPLAY_MATH_BLOCK.finditer(text):
            prefix = text[cursor:match.start()].strip()
            if prefix:
                SSAS._publish(Markdown(prefix))
            math_body = match.group(1).strip()
            if math_body:
                SSAS._publish(_DisplayMath(math_body))
                rendered_math = True
            cursor = match.end()
        suffix = text[cursor:].strip()
        if suffix:
            SSAS._publish(Markdown(suffix))
        elif not rendered_math:
            SSAS._publish(Markdown(text))

    @staticmethod
    def _markdown_type():
        try:
            from IPython.display import Markdown
        except ImportError:
            return None
        return Markdown

    @staticmethod
    def _publish(obj) -> None:
        try:
            from IPython.core.interactiveshell import InteractiveShell
        except ImportError:
            return False
        if not InteractiveShell.initialized():
            return
        shell = InteractiveShell.instance()
        data = {}
        markdown = getattr(obj, "_repr_markdown_", None)
        latex = getattr(obj, "_repr_latex_", None)
        if markdown is not None:
            data["text/markdown"] = markdown()
        if latex is not None:
            data["text/latex"] = latex()
        if not data:
            data["text/plain"] = repr(obj)
        shell.display_pub.publish(data, metadata={})

    @staticmethod
    def _safe_stem(text: str) -> str:
        cleaned = re.sub(r"[^0-9A-Za-z]+", "_", text).strip("_").lower()
        return cleaned[:80] or "response"

    def _route_path(self, path: str | PathLike[str]) -> Path:
        candidate = Path(path).expanduser()
        if self.absolute_route:
            return candidate.resolve()
        if candidate.is_absolute():
            return candidate
        return self.route_dir / candidate

    @staticmethod
    def _caller_dir() -> Path:
        this_file = Path(__file__).resolve()
        cwd = Path.cwd().resolve()
        frames = inspect.stack()[2:]
        if any(frame.filename and frame.filename.startswith("<") for frame in frames):
            return cwd
        external_candidates: list[Path] = []
        for frame in frames:
            filename = frame.filename
            if not filename or filename.startswith("<"):
                continue
            path = Path(filename).resolve()
            if path == this_file or not path.exists() or not path.is_file():
                continue
            if path.is_relative_to(cwd):
                return path.parent
            if "site-packages" not in path.parts and "dist-packages" not in path.parts:
                external_candidates.append(path)
        if external_candidates:
            return external_candidates[0].parent
        return cwd
