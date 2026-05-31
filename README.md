# small-signal-equation-extractor

`small-signal-equation-extractor` extracts symbolic small-signal equations from a restricted SPICE-like netlist. It builds a linear modified nodal analysis system, solves exact rational responses in the Laplace variable `s`, and writes the result as Markdown equations.

The project is an equation extractor, not a circuit simulator. It does not compute an operating point, parse process design kit model cards, or evaluate Berkeley Short-channel Insulated-Gate Field-Effect Transistor Model version 4 (BSIM4) code. Compact-device coefficients are symbolic unless a netlist element explicitly provides a numeric controlled-source gain.

## Installation

Python 3.10 or later is required.

Clone the repository and run it from the repository root:

```bash
git clone https://github.com/YoungchanJo/small-signal-equation-extractor.git
cd small-signal-equation-extractor
```

The core extractor uses only the Python standard library. IPython is optional and is used only for notebook display with `show=True`.

## Quick start

Use `SSAS.solve(...)`. A request is a voltage probe, a current probe, or one ratio of two probes. Every successful call writes one Markdown file and returns `None`.

```python
from SSAS import SSAS

ssas = SSAS(output_dir="tutorial_outputs")

netlist = """
V1 in 0 ac=1
R1 in out 1
C1 out 0 1
"""

ssas.solve(
    netlist,
    "V(out) / V(in)",
    "rc_lowpass.md",
    source="inline",
    show=True,
)
```

$$H(s)=\frac{1}{s R_{1} C_{1}+1}$$

For a file-based netlist, use `source="file"`:

```python
ssas.solve(
    "netlists/cs_amp.sp",
    "V(out) / V(in)",
    "cs_gain.md",
    source="file",
    eq_circuit=True,
    body="on",
    show=True,
)
```

Relative `output_dir`, output filenames, and file netlist paths are resolved from the caller `.py` file directory. In notebooks, they are resolved from the notebook working directory, which should be the directory containing `tutorial.ipynb`. Use `absolute_route=True` only when absolute-path routing is intended.

If the same output filename is used again, the file is overwritten.

## Request syntax

```text
V(node)
I(device)
V(node1) / V(node2)
V(node) / I(device)
I(device1) / V(node)
I(device1) / I(device2)
```

`I(device)` is available only for devices with a public two-terminal current convention. It is not a public current probe for MOSFETs, bipolar junction transistors, controlled sources, or subcircuit instances.

## Compact-device stamp modes

`eq_circuit=True` selects equivalent-circuit mode for supported compact devices. MOSFETs use equivalent-circuit symbols such as `g_m`, `g_mb`, `r_o`, and terminal-pair capacitances. Bipolar junction transistors use their equivalent-circuit transconductance, input resistance, output resistance, and terminal-pair capacitances.

`eq_circuit=False` selects coefficient mode. For MOSFETs, this is a BSIM4-compatible external-port coefficient stamp: a separate wrapper may evaluate the device model and supply terminal-current and terminal-charge derivative coefficients after reducing internal model nodes to the external port. For bipolar junction transistors, coefficient mode uses symbolic differential conductance and device differential capacitance coefficients.

In both modes, the extractor stamps the symbolic coefficients it is given. It does not derive compact-model coefficients from nonlinear device equations.

## Main options

| Option | Meaning |
|---|---|
| `source` | Required. Use `"inline"` for netlist text and `"file"` for a netlist path. |
| `eq_circuit` | Selects equivalent-circuit mode when `True` and coefficient mode when `False`. |
| `body` | Enables MOSFET body terms when set to `"on"`. |
| `zero_pole` | Appends analytic zero-pole post-analysis after the ordinary response expression. |
| `fraction_symbol_threshold` | Controls when long coefficient expressions are replaced by auxiliary `E_i` definitions. |
| `off_device` | Removes selected compact-device conductance or current terms. |
| `off_device_cap` | Removes selected compact-device capacitance terms. |
| `zero_cap` | Removes selected passive, equivalent-circuit, or coefficient-mode capacitance terms. |
| `zero_resistance` | Shorts selected passive or equivalent-circuit resistances. |
| `zero_inductance` | Shorts selected inductances. |
| `equal_passive` | Aliases equal passive resistor, capacitor, or inductor values only in the final expression. |
| `subckt` | Selects one `.subckt` definition for extraction. |
| `subckt_boundary` | Controls skipped subcircuit boundary treatment. |
| `show` | Displays generated Markdown in an IPython notebook. |

Ideal switch primitives are not supported; SPICE switch-style `S` and `W` elements are outside the implemented netlist subset.

## Current implementation limits

The present MOSFET coefficient-mode implementation supports only quasi-static coefficient stamps. Here, quasi-static means that the device contribution is represented by terminal-current coefficients and terminal-charge derivative coefficients, so the Laplace-domain contribution is affine in `s`.

The extractor does not implement transient non-quasi-static MOSFET states, frequency-domain non-quasi-static effective admittance, BSIM4 noise-source mapping, or ideal switch primitives. It should be used as a quasi-static symbolic external-port coefficient extractor until those interfaces are added.

## Documentation

- `USAGE.md`: API usage, output-oriented option examples, and netlist subset.
- `mna_device_equations.md`: modified nodal analysis stamps and sign conventions.
- `tutorial.ipynb`: executable tutorial with rendered outputs.
