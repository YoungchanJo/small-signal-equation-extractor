# Usage

`SSAS` builds symbolic small-signal responses from a restricted SPICE-like netlist. It parses the netlist, builds a linear modified nodal analysis system, solves one requested response, and saves one Markdown file. Every successful `solve(...)` call returns `None`.

The extractor is not a simulator. It does not compute an operating point, parse process design kit model cards, or evaluate Berkeley Short-channel Insulated-Gate Field-Effect Transistor Model version 4 (BSIM4) code. Compact-device coefficients are symbolic unless a netlist element explicitly provides a numeric controlled-source gain.

## 1. Basic workflow

Use one public method: `solve(...)`.

Minimal calls are:

- Inline netlist: `ssas.solve(netlist, "V(out) / V(in)", "rc_lowpass.md", source="inline")`
- File netlist: `ssas.solve("netlists/cs_amp.sp", "V(out) / V(in)", "cs_gain.md", source="file")`

`source` is mandatory. Use `source="inline"` for netlist text and `source="file"` for a netlist path. `SSAS` does not infer whether the first argument is inline netlist text or a path.

Representative output fragment:

$$
H(s)=\frac{1}{C_{1}R_{1}s+1}
$$

## 2. Path routing

Construct `SSAS` with `SSAS(output_dir="tutorial_outputs")` for ordinary relative routing.

With the default `absolute_route=False`, relative paths are routed from the caller location: in a `.py` script, from the directory containing that script; in a notebook, from the notebook working directory, which should be the directory containing `tutorial.ipynb`.

The same route rule applies to `output_dir`, relative output filenames passed to `solve(...)`, and relative netlist paths passed with `source="file"`.

Use `absolute_route=True` only when the paths passed to `SSAS` and `solve(...)` are intended to be absolute-path based, for example `SSAS(output_dir="/tmp/ssas_outputs", absolute_route=True)`.

If an output file already exists, the new result overwrites it.

## 3. Request syntax

A request is one probe or one top-level ratio of two probes: `V(node)`, `I(device)`, `V(node1) / V(node2)`, `V(node) / I(device)`, `I(device1) / V(node)`, or `I(device1) / I(device2)`.

For example, `"V(out)"` solves a single node voltage, `"I(RD)"` solves a passive-device current, and `"I(RD) / V(in)"` solves a transconductance-like ratio.

Representative output changes:

$$
\begin{aligned}
\text{voltage probe:}\quad & H(s)=V_{out}\text{ expression},\\
\text{current probe:}\quad & H(s)=I_{RD}\text{ expression}.
\end{aligned}
$$

The left-hand side is always `H(s)`; the request controls the solved expression on the right-hand side. `V(0)` and `V(gnd)` are valid zero-voltage probes. MOSFETs, bipolar junction transistors, controlled sources, and subcircuit instances are not public current probes.

## 4. Rendering options

### `zero_pole`

Use `zero_pole=True` to append zero-pole analysis after the ordinary response. It does not replace `H(s)`.

Without `zero_pole`, the file ends after the response and its auxiliary definitions:

$$
H(s)=\frac{\cdots}{\cdots}
$$

With `zero_pole=True`, the same response is kept and a final section is appended:

$$
H(s)=\frac{\cdots}{\cdots}
$$

$$
\text{zero-pole analysis}
$$

$$
Z_{1}=s-\left(\cdots\right),
\qquad
P_{1}=s-\left(\cdots\right)
$$

If a zero or pole factor cannot be solved analytically, the Markdown reports the status instead of raising an exception:

$$
\text{remaining poles: analytic extraction failed}
$$

### `fraction_symbol_threshold`

The renderer always simplifies the displayed rational expression when simplification reduces rendering cost. For example, a conductance-first internal form may be displayed as a resistance-based compact expression when that form is shorter.

`fraction_symbol_threshold` controls only when a long coefficient expression is replaced by an auxiliary definition such as `E_{1}`. The replaced expression can contain products, sums, powers, or rational subexpressions. This option does not enable or disable ordinary expression simplification. The default is `10`. Use `fraction_symbol_threshold=0` to suppress `E_i` coefficient replacement.

Typical output when no coefficient replacement is used:

$$
H(s)=\frac{a_{2}s^{2}+a_{1}s+a_{0}}{b_{2}s^{2}+b_{1}s+b_{0}}
$$

With a small positive threshold, long coefficient expressions are lifted into definitions, but the transfer function is otherwise the same response:

$$
H(s)=\frac{E_{1}s+E_{2}}{E_{3}s^{2}+E_{4}s+E_{5}}
$$

$$
E_{1}=\cdots
$$

### `show`

Use `show=True` in a notebook to display the generated Markdown immediately. It does not change the saved file content.

## 5. Circuit-mode options

### `eq_circuit`

`eq_circuit=True` selects equivalent-circuit mode for supported compact devices. MOSFETs use `g_m`, optional `g_mb`, `r_o`, and terminal-pair capacitances such as `C_{mn1,gs}` and `C_{mn1,gd}`. Bipolar junction transistors use an equivalent-circuit stamp with transconductance, input resistance, output resistance, and terminal-pair capacitances.

Representative MOSFET equivalent-circuit fragment:

$$
H(s)=\frac{-g_{mn1}r_{mn1}R_{D}}{r_{mn1}+R_{D}+\cdots}
$$

`eq_circuit=False` selects coefficient mode for supported compact devices. For MOSFETs, this is a BSIM4-compatible external-port coefficient stamp: a separate BSIM4 wrapper may supply terminal-current and terminal-charge derivative coefficients after evaluating BSIM4 and reducing internal BSIM4 nodes to the external terminal port. For bipolar junction transistors, this mode uses symbolic differential conductance and device differential capacitance coefficients. The extractor itself only stamps the supplied symbolic coefficients.

Representative MOSFET coefficient-mode fragment:

$$
H(s)=\frac{-g_{mn1,gs,d}+sC_{mn1,g,d}+\cdots}{g_{mn1,ds,d}+sC_{mn1,d,d}+\cdots}
$$

Representative bipolar junction transistor mode distinction:

$$
\begin{aligned}
\text{equivalent-circuit mode:}\quad & g_{q1},\ r_{q1,o},\ r_{q1,\pi},\ C_{q1,bc}\text{ may appear},\\
\text{coefficient mode:}\quad & g_{q1,be,c},\ g_{q1,ce,c},\ C_{q1,b,c}\text{ may appear}.
\end{aligned}
$$

### `body`

Use `body="off"` to omit MOSFET body rows and body-dependent terms from the MOSFET stamp. Use `body="on"` to include body-dependent conductance and device differential capacitance terms.

Representative change in equivalent-circuit mode:

$$
\begin{aligned}
\text{body terms omitted:}\quad & H(s)=\frac{\cdots g_{mn1}\cdots}{\cdots},\\
\text{body terms included:}\quad & H(s)=\frac{\cdots g_{mn1}+\cdots g_{mb,mn1}\cdots}{\cdots}.
\end{aligned}
$$

Representative change in MOSFET coefficient mode:

$$
\begin{aligned}
\text{body terms omitted:}\quad & \text{body-dependent coefficient rows are omitted},\\
\text{body terms included:}\quad & \text{terms such as }g_{mn1,bs,d}\text{ and }C_{mn1,b,d}\text{ may appear}.
\end{aligned}
$$

The option is meaningful for MOSFETs in both `eq_circuit=True` and `eq_circuit=False` modes. It does not affect bipolar junction transistor stamps.

## 6. Term-selection options

The following options apply to compact-device instances in both circuit modes unless stated otherwise. The examples below use MOSFET and bipolar junction transistor symbols.

### `off_device`

`off_device` removes conductance or current terms for selected compact-device instances. It leaves device capacitance terms enabled. Use `off_device="mn1"` or `off_device=["mn1", "mp3"]`.

Representative equivalent-circuit change:

$$
\begin{aligned}
\text{baseline:}\quad & H(s)=\frac{\cdots g_{mn1}\cdots r_{mn1}\cdots C_{mn1,gd}\cdots}{\cdots},\\
\text{selected conductance terms removed:}\quad & H(s)=\frac{\cdots C_{mn1,gd}\cdots}{\cdots}.
\end{aligned}
$$

Representative coefficient-mode change:

$$
\begin{aligned}
\text{baseline:}\quad & g_{mn1,gs,d},\ g_{mn1,ds,d},\ C_{mn1,g,d}\text{ may appear},\\
\text{selected conductance terms removed:}\quad & g_{mn1,\ast,\ast}\text{ terms are removed; }C_{mn1,\ast,\ast}\text{ terms remain}.
\end{aligned}
$$

### `off_device_cap`

`off_device_cap` removes capacitance terms for selected compact-device instances. It leaves conductance or current terms enabled. Use `off_device_cap="mn1"` or `off_device_cap=["mn1", "q2"]`.

Representative equivalent-circuit change:

$$
\begin{aligned}
\text{baseline:}\quad & g_{mn1},\ r_{mn1},\ C_{mn1,gd}\text{ may appear},\\
\text{selected capacitance terms removed:}\quad & g_{mn1},\ r_{mn1}\text{ remain; }C_{mn1,\ast}\text{ terms are removed}.
\end{aligned}
$$

Representative coefficient-mode change:

$$
\begin{aligned}
\text{baseline:}\quad & g_{mn1,gs,d},\ C_{mn1,g,d}\text{ may appear},\\
\text{selected capacitance terms removed:}\quad & g_{mn1,gs,d}\text{ remains; }C_{mn1,\ast,\ast}\text{ terms are removed}.
\end{aligned}
$$

`off_device` and `off_device_cap` are independent:

| Options | Conductance/current terms | Capacitance terms |
|---|---|---|
| neither option | included | included |
| `off_device="mn1"` | removed | included |
| `off_device_cap="mn1"` | included | removed |
| both contain the same device | removed | removed |

### `zero_cap`

`zero_cap` removes selected capacitance terms by symbol name. Use `zero_cap="C_mn1_gd"` for an equivalent-circuit terminal-pair capacitance, `zero_cap="C_mn1_g_d"` for a MOSFET coefficient-mode device differential capacitance coefficient, `zero_cap="C_q1_b_c"` for a bipolar junction transistor coefficient-mode device differential capacitance coefficient, or `zero_cap="C_1"` for a passive capacitor.

Equivalent-circuit MOSFET example:

$$
\begin{aligned}
\text{baseline:}\quad & C_{mn1,gd},\ C_{mn1,ds},\ C_{mn1,db}\text{ may appear},\\
\text{selected terminal-pair capacitance removed:}\quad & C_{mn1,gd}\text{ is removed; }C_{mn1,ds}\text{ and }C_{mn1,db}\text{ remain}.
\end{aligned}
$$

MOSFET coefficient-mode example:

$$
\begin{aligned}
\text{baseline:}\quad & C_{mn1,g,d},\ C_{mn1,d,d}\text{ may appear},\\
\text{selected differential capacitance removed:}\quad & C_{mn1,g,d}\text{ is removed; }C_{mn1,d,d}\text{ remains}.
\end{aligned}
$$

Passive capacitor example:

$$
\begin{aligned}
\text{baseline:}\quad & H(s)=\frac{1}{C_{1}R_{1}s+1},\\
\text{selected passive capacitance removed:}\quad & H(s)=1.
\end{aligned}
$$

### `zero_resistance`

`zero_resistance` replaces selected passive or equivalent-circuit resistance symbols with ideal shorts. Use `zero_resistance="R_D"` for a passive resistor or `zero_resistance="r_mn1"` for an equivalent-circuit compact-device output resistance.

Representative change for a common-source load resistor:

$$
\begin{aligned}
\text{baseline:}\quad & H(s)=\frac{\cdots R_{D}\cdots}{\cdots},\\
\text{selected resistance shorted:}\quad & H(s)=0.
\end{aligned}
$$

Equivalent MOSFET or bipolar junction transistor internal resistance symbols exist only in equivalent-circuit mode. Passive-resistance shorting applies in both `eq_circuit=True` and `eq_circuit=False`.

### `zero_inductance`

`zero_inductance` replaces selected inductance symbols with ideal shorts. Use `zero_inductance="L_1"` for a passive inductor.

Representative change:

$$
\begin{aligned}
\text{baseline:}\quad & H(s)=\frac{L_{1}s}{L_{1}s+R_{1}},\\
\text{selected inductance shorted:}\quad & H(s)=0.
\end{aligned}
$$

This option is independent of `eq_circuit`; it applies to passive inductors.

### `equal_passive`

`equal_passive` aliases passive resistor, capacitor, or inductor value symbols only in the final rendered expression. It does not change topology, stamping, zero-element handling, or current-probe names. Use `equal_passive=("R_X", ["r1", "r2"])` for one group or `equal_passive=[("R_X", ["r1", "r2"]), ("C_X", ["c1", "c2"])]` for multiple groups.

Representative change:

$$
\begin{aligned}
\text{baseline:}\quad & H(s)=\frac{\cdots R_{1}+\cdots R_{2}\cdots}{\cdots},\\
\text{selected passive values aliased:}\quad & H(s)=\frac{\cdots R_{X}+\cdots R_{X}\cdots}{\cdots}.
\end{aligned}
$$

Each group must contain only one passive type. MOSFETs, bipolar junction transistors, sources, and mixed resistor/capacitor/inductor groups are not valid `equal_passive` targets.

## 7. Subcircuits

Top-level `X...` instances are ignored unless a `subckt` name is selected. When `subckt` is selected, the first matching top-level instance is mapped into the main circuit. Use `subckt="AMP"` to select a subcircuit definition.

Representative output effect:

$$
\begin{aligned}
\text{no selected subcircuit:}\quad & \text{top-level selected circuit is solved},\\
\text{selected subcircuit expanded:}\quad & \text{the selected subcircuit instance is expanded into the solved circuit}.
\end{aligned}
$$

`subckt_boundary="ground"` ties boundary nodes of skipped sibling subcircuits to small-signal ground. `subckt_boundary="symbolic"` exposes them as independent symbolic voltages.

Representative output effect:

$$
\begin{aligned}
\text{grounded skipped-boundary nodes:}\quad & \text{skipped boundary voltages are treated as }0,\\
\text{symbolic skipped-boundary nodes:}\quad & \text{skipped boundary voltages may appear as symbolic inputs}.
\end{aligned}
$$

## 8. Netlist subset

Supported primitives:

| Prefix | Device |
|---|---|
| `V` | independent voltage source |
| `I` | independent current source |
| `R` | resistor |
| `C` | capacitor |
| `L` | inductor |
| `D` | diode |
| `E` | voltage-controlled voltage source |
| `G` | voltage-controlled current source |
| `F` | current-controlled current source |
| `H` | current-controlled voltage source |
| `M`, `MN`, `MP`, `NM`, `PM` | MOSFET symbolic small-signal stamp |
| `Q`, `QN`, `QP`, `NQ`, `PQ` | bipolar junction transistor symbolic small-signal stamp |

Ideal switches are not supported. SPICE switch-style primitives such as `S` or `W` are outside the implemented netlist subset.

Unsupported SPICE features are ignored or rejected. The extractor does not implement `.model`, `.param`, `.include`, nonlinear analysis directives, measurements, or simulator-specific control blocks.
