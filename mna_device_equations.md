# Modified Nodal Analysis Device Equations Used by the Extractor

This document describes the equations used when `SSAS` builds the symbolic modified nodal analysis (MNA) system. It is a stamp reference for the extractor, not a compact-model derivation. Metal-oxide-semiconductor field-effect transistor (MOSFET), bipolar junction transistor (BJT), and diode compact-model coefficients are treated as symbolic small-signal parameters.

Minimal call pattern:

```python
from SSAS import SSAS

ssas = SSAS(output_dir="tutorial_outputs")
ssas.solve(
    "netlists/cs_amp.sp",
    "V(out) / V(in)",
    "cs_gain.md",
    source="file",
    eq_circuit=True,
    body="on",
)
```

The `eq_circuit`, `body`, `off_device`, `off_device_cap`, `zero_cap`, `zero_resistance`, and `zero_inductance` options select which stamp terms below are included.

## 1. Conventions

### 1.1 MNA System

Each row is stored as a residual equation:

$$
F_k(x,s)=0
$$

The assembled linear system has the form:

$$
A(s)x+c(s)=0
$$

where `x` contains unknown node voltages and branch currents. Symbolic source voltages, source currents, device parameters, and the Laplace variable `s` appear in the coefficients or constant terms.

### 1.2 Node Voltages

Node voltage is written as:

$$
v_n
$$

The voltage from node `a` to node `b` is:

$$
v_{ab}=v_a-v_b
$$

Ground node names are `0` and `gnd`, and their voltage is zero:

$$
v_0=0
$$

Ground boundary nodes also have zero voltage. Symbolic boundary nodes are skipped by KCL and exposed as independent symbols `v_{node}`.

### 1.3 KCL Sign Convention

A device current from node `a` to node `b` is written as:

$$
i_{ab}:a\rightarrow b
$$

The extractor writes KCL rows as sums of currents leaving each node:

$$
F_a \leftarrow F_a+i_{ab}
$$

$$
F_b \leftarrow F_b-i_{ab}
$$

Thus every terminal-current expression below is the current contribution added to the KCL row of that terminal. KCL rows are not created for ground nodes, ground boundary nodes, or symbolic boundary nodes.

### 1.4 Laplace Variable

Capacitance and inductance stamps are written in the Laplace domain:

$$
\frac{d}{dt}\mapsto s
$$

For example:

$$
i_C=sC(v_a-v_b)
$$

### 1.5 Symbol Names

The equations use mathematical names such as `R_D` and `C_{mn1,gd}`. In Python options and numeric assignments, those same names are passed through plain normalization: braces and punctuation become underscores. For example, `C_{mn1,gd}` is addressed as `C_mn1_gd`, and `r_{q1,o}` is addressed as `r_q1_o`.

Small-signal element values are symbolic unless a controlled-source gain is written as a numeric token. Passive resistor, capacitor, and inductor value tokens are required by syntax, but they name symbolic parameters rather than substituting numeric values.

## 2. Independent Sources

### 2.1 Voltage Source `V`

Netlist form:

```spice
Vname n+ n- ...
```

The voltage source introduces a branch-current unknown:

$$
i_V:n_+\rightarrow n_-
$$

KCL stamp:

$$
F_{n_+}\leftarrow F_{n_+}+i_V
$$

$$
F_{n_-}\leftarrow F_{n_-}-i_V
$$

If the AC value is missing or zero, the source is a small-signal short:

$$
v_{n_+}-v_{n_-}=0
$$

If the AC value is nonzero and the source is ground-referenced, the non-ground node is exposed as an independent symbolic voltage `v_{node}`. No voltage-source constraint is added in that case. The AC value is used only to detect nonzero excitation, not as a numeric amplitude.

If the AC value is nonzero and the source is not ground-referenced, the constraint is:

$$
v_{n_+}-v_{n_-}-V_V=0
$$

where `V_V` is a symbolic input voltage derived from the source instance name.

### 2.2 Current Source `I`

Netlist form:

```spice
Iname n+ n- ...
```

The source current direction is:

$$
i_I:n_+\rightarrow n_-
$$

If the AC value is missing or zero, no current is stamped and the registered current probe value is zero:

$$
i_I=0
$$

If the AC value is nonzero, the source stamps a symbolic input current:

$$
F_{n_+}\leftarrow F_{n_+}+i_I
$$

$$
F_{n_-}\leftarrow F_{n_-}-i_I
$$

As with voltage sources, the AC value only distinguishes zero from nonzero.

## 3. Controlled Sources

Controlled-source gain tokens may be numeric or symbolic. Numeric gains are stored as exact rational constants when possible. Symbolic gains become ordinary symbolic coefficients.

### 3.1 VCVS `E`

Netlist form:

```spice
Ename n+ n- nc+ nc- gain
```

An output branch-current unknown is added:

$$
i_E:n_+\rightarrow n_-
$$

KCL stamp:

$$
F_{n_+}\leftarrow F_{n_+}+i_E
$$

$$
F_{n_-}\leftarrow F_{n_-}-i_E
$$

Constraint:

$$
v_{n_+}-v_{n_-}-A_E(v_{nc+}-v_{nc-})=0
$$

### 3.2 VCCS `G`

Netlist form:

```spice
Gname n+ n- nc+ nc- gain
```

Current direction:

$$
i_G:n_+\rightarrow n_-
$$

Current expression:

$$
i_G=A_G(v_{nc+}-v_{nc-})
$$

KCL stamp:

$$
F_{n_+}\leftarrow F_{n_+}+A_G(v_{nc+}-v_{nc-})
$$

$$
F_{n_-}\leftarrow F_{n_-}-A_G(v_{nc+}-v_{nc-})
$$

### 3.3 CCCS `F`

Netlist form:

```spice
Fname n+ n- ctrl_source gain
```

Current direction:

$$
i_F:n_+\rightarrow n_-
$$

The control quantity is the branch current associated with `ctrl_source`:

$$
i_F=A_F i_{\mathrm{ctrl}}
$$

KCL stamp:

$$
F_{n_+}\leftarrow F_{n_+}+A_F i_{\mathrm{ctrl}}
$$

$$
F_{n_-}\leftarrow F_{n_-}-A_F i_{\mathrm{ctrl}}
$$

`ctrl_source` should name a source/device with a meaningful branch-current unknown. `SSAS` does not elaborate arbitrary subcircuits to create such a current.

### 3.4 CCVS `H`

Netlist form:

```spice
Hname n+ n- ctrl_source gain
```

An output branch-current unknown is added:

$$
i_H:n_+\rightarrow n_-
$$

KCL stamp:

$$
F_{n_+}\leftarrow F_{n_+}+i_H
$$

$$
F_{n_-}\leftarrow F_{n_-}-i_H
$$

Constraint:

$$
v_{n_+}-v_{n_-}-A_H i_{\mathrm{ctrl}}=0
$$

## 4. Passive Resistor, Capacitor, and Inductor Stamps

R, C, and L values are required syntactically but are not substituted numerically. The symbolic parameter name is derived from the device name. For example, `RD d 0 10k` uses `R_D`, and `Cload out 0 1p` uses `C_{load}`.

### 4.1 Resistor `R`

Netlist form:

```spice
Rname n+ n- value
```

Current direction:

$$
i_R:n_+\rightarrow n_-
$$

Current expression:

$$
i_R=\frac{v_{n_+}-v_{n_-}}{R_R}
$$

KCL stamp:

$$
F_{n_+}\leftarrow F_{n_+}+\frac{v_{n_+}-v_{n_-}}{R_R}
$$

$$
F_{n_-}\leftarrow F_{n_-}-\frac{v_{n_+}-v_{n_-}}{R_R}
$$

### 4.2 Capacitor `C`

Netlist form:

```spice
Cname n+ n- value
```

Current direction:

$$
i_C:n_+\rightarrow n_-
$$

Current expression:

$$
i_C=sC_C(v_{n_+}-v_{n_-})
$$

KCL stamp:

$$
F_{n_+}\leftarrow F_{n_+}+sC_C(v_{n_+}-v_{n_-})
$$

$$
F_{n_-}\leftarrow F_{n_-}-sC_C(v_{n_+}-v_{n_-})
$$

### 4.3 Inductor `L`

Netlist form:

```spice
Lname n+ n- value
```

The inductor introduces a branch-current unknown:

$$
i_L:n_+\rightarrow n_-
$$

KCL stamp:

$$
F_{n_+}\leftarrow F_{n_+}+i_L
$$

$$
F_{n_-}\leftarrow F_{n_-}-i_L
$$

Constraint:

$$
v_{n_+}-v_{n_-}-sL_L i_L=0
$$

## 5. Diode Stamp

Netlist form:

```spice
Dname n+ n- ...
```

The first terminal defines the positive current direction. The model name and remaining diode parameters are accepted syntactically but are not evaluated by the extractor.

For canonical diode name `D`, the symbolic small-signal current is:

$$
i_D=(g_D+sC_D)(v_{n+}-v_{n-}),\qquad i_D:n_+\rightarrow n_-
$$

KCL stamp:

$$
F_{n+}\leftarrow F_{n+}+i_D
$$

$$
F_{n-}\leftarrow F_{n-}-i_D
$$

The diode stamp is independent of `eq_circuit`.

## 6. MOSFET Stamps

MOSFET netlist terminal order is:

```spice
Mname drain gate source body ...
```

The equations use local terminal labels:

$$
d=\mathrm{drain},\quad g=\mathrm{gate},\quad s=\mathrm{source},\quad b=\mathrm{body}
$$

`eq_circuit=False` selects coefficient mode for compact devices. For MOSFETs, this is a BSIM4-compatible external-port coefficient stamp. `eq_circuit=True` selects equivalent-circuit mode for MOSFETs and bipolar junction transistors. Diodes use the same two-terminal symbolic stamp in both modes.

The modes use different meanings for MOSFET `C` symbols. In equivalent-circuit mode, a `C` symbol is a two-terminal Meyer-style capacitor. In MOSFET coefficient mode, terminal charges are evaluated first; capacitance coefficients are charge derivatives, not standalone circuit capacitors.

### 6.1 MOSFET coefficient stamp: `eq_circuit=False`

This mode stamps an external-port coefficient relation. A BSIM4 wrapper may supply these coefficients after evaluating BSIM4, normalizing external terminal currents and charges, and reducing BSIM4-internal nodes. The extractor itself only stamps the supplied symbolic coefficients; it does not evaluate BSIM4.

The active terminal row set is:

$$
\mathcal{R}_M=\{g,d,s\}
$$

when `body="off"`, and:

$$
\mathcal{R}_M=\{g,d,s,b\}
$$

when `body="on"`.

For each active row terminal `r`, the extractor adds:

$$
F_r\leftarrow F_r+i_{M,r}
$$

if the circuit node connected to that terminal has a KCL row.

Conductance part:

$$
i_{M,r}^{(G)}=
g_{M,gs,r}(v_g-v_s)+g_{M,ds,r}(v_d-v_s)
$$

With `body="on"`:

$$
i_{M,r}^{(G)}=
g_{M,gs,r}(v_g-v_s)+g_{M,ds,r}(v_d-v_s)+g_{M,bs,r}(v_b-v_s)
$$

Charge part:

$$
i_{M,r}^{(Q)}=s\sum_{c\in\mathcal{R}_M} C_{M,c,r}v_c
$$

Total current contribution:

$$
i_{M,r}=i_{M,r}^{(G)}+i_{M,r}^{(Q)}
$$

The final suffix in `g_{M,gs,r}` or `C_{M,c,r}` is the KCL row terminal. The `C_{M,c,r}` symbols are terminal-charge Jacobian entries. They are not terminal-pair capacitors.

For `body="on"`, the compact matrix form is:

$$
\boldsymbol{i}_{M}=G_M\boldsymbol{u}_{M}+sC_M\boldsymbol{v}_{M}
$$

where the vectors are defined as follows:

$$
\boldsymbol{i}_{M}=\left(i_{M,g},i_{M,d},i_{M,s},i_{M,b}\right)^{T}
$$

$$
\boldsymbol{u}_{M}=\left(v_g-v_s,\ v_d-v_s,\ v_b-v_s\right)^{T}
$$

$$
\boldsymbol{v}_{M}=\left(v_g,v_d,v_s,v_b\right)^{T}
$$

The rows of $G_M$ are ordered as $g,d,s,b$, and its columns are ordered as $gs,ds,bs$:

| row | $gs$ | $ds$ | $bs$ |
|---|---|---|---|
| $g$ | $g_{M,gs,g}$ | $g_{M,ds,g}$ | $g_{M,bs,g}$ |
| $d$ | $g_{M,gs,d}$ | $g_{M,ds,d}$ | $g_{M,bs,d}$ |
| $s$ | $g_{M,gs,s}$ | $g_{M,ds,s}$ | $g_{M,bs,s}$ |
| $b$ | $g_{M,gs,b}$ | $g_{M,ds,b}$ | $g_{M,bs,b}$ |

The rows of $C_M$ are ordered as $g,d,s,b$, and its columns are ordered as $g,d,s,b$:

| row | $g$ | $d$ | $s$ | $b$ |
|---|---|---|---|---|
| $g$ | $C_{M,g,g}$ | $C_{M,d,g}$ | $C_{M,s,g}$ | $C_{M,b,g}$ |
| $d$ | $C_{M,g,d}$ | $C_{M,d,d}$ | $C_{M,s,d}$ | $C_{M,b,d}$ |
| $s$ | $C_{M,g,s}$ | $C_{M,d,s}$ | $C_{M,s,s}$ | $C_{M,b,s}$ |
| $b$ | $C_{M,g,b}$ | $C_{M,d,b}$ | $C_{M,s,b}$ | $C_{M,b,b}$ |

For `body="off"`, remove the `b` row, `b` charge column, and `bs` conductance column.

### 6.2 MOSFET equivalent-circuit stamp: `eq_circuit=True`

This mode stamps a MOSFET equivalent-circuit relation. It contains:

- a controlled current source from drain to source,
- an optional body transconductance source from drain to source,
- an output resistance between drain and source,
- terminal-pair capacitors.

These capacitors are Meyer-style equivalent capacitors. They should not be identified with MOSFET coefficient-mode `C` terms or BSIM4 charge derivatives.

Conductance/current terms:

$$
i_{gm}=g_M(v_g-v_s),\qquad i_{gm}:d\rightarrow s
$$

$$
i_{ro}=\frac{v_d-v_s}{r_M},\qquad i_{ro}:d\rightarrow s
$$

When `body="on"`:

$$
i_{gmb}=g_{M,mb}(v_b-v_s),\qquad i_{gmb}:d\rightarrow s
$$

Capacitance pairs for `body="off"`:

$$
(g,s),\quad (g,d),\quad (d,s)
$$

Additional pairs for `body="on"`:

$$
(g,b),\quad (d,b),\quad (s,b)
$$

Each capacitance current is directed from the first terminal in the pair to the second:

$$
i_{Cxy}=sC_{M,xy}(v_x-v_y)
$$

Terminal current contributions for `body="off"`:

$$
i_g=sC_{M,gs}(v_g-v_s)+sC_{M,gd}(v_g-v_d)
$$

$$
i_d=g_M(v_g-v_s)+\frac{v_d-v_s}{r_M}
+sC_{M,gd}(v_d-v_g)+sC_{M,ds}(v_d-v_s)
$$

$$
i_s=-g_M(v_g-v_s)-\frac{v_d-v_s}{r_M}
+sC_{M,gs}(v_s-v_g)+sC_{M,ds}(v_s-v_d)
$$

Terminal current contributions for `body="on"`:

$$
i_g=sC_{M,gs}(v_g-v_s)+sC_{M,gd}(v_g-v_d)+sC_{M,gb}(v_g-v_b)
$$

$$
i_d=g_M(v_g-v_s)+g_{M,mb}(v_b-v_s)+\frac{v_d-v_s}{r_M}
+sC_{M,gd}(v_d-v_g)+sC_{M,ds}(v_d-v_s)+sC_{M,db}(v_d-v_b)
$$

$$
i_s=-g_M(v_g-v_s)-g_{M,mb}(v_b-v_s)-\frac{v_d-v_s}{r_M}
+sC_{M,gs}(v_s-v_g)+sC_{M,ds}(v_s-v_d)+sC_{M,sb}(v_s-v_b)
$$

$$
i_b=sC_{M,gb}(v_b-v_g)+sC_{M,db}(v_b-v_d)+sC_{M,sb}(v_b-v_s)
$$

## 7. BJT Stamps

BJT netlist terminal order is:

```spice
Qname collector base emitter ...
```

The equations use local terminal labels:

$$
c=\mathrm{collector},\quad b=\mathrm{base},\quad e=\mathrm{emitter}
$$

### 7.1 Compact BJT: `eq_circuit=False`

The terminal row set is:

$$
\mathcal{R}_Q=\{c,b,e\}
$$

For each row terminal `r`, the extractor adds:

$$
F_r\leftarrow F_r+i_{Q,r}
$$

Conductance part:

$$
i_{Q,r}^{(G)}=
g_{Q,be,r}(v_b-v_e)+g_{Q,ce,r}(v_c-v_e)
$$

Capacitance part:

$$
i_{Q,r}^{(C)}=s\sum_{k\in\{c,b,e\}} C_{Q,k,r}v_k
$$

Total compact current:

$$
i_{Q,r}=
g_{Q,be,r}(v_b-v_e)+g_{Q,ce,r}(v_c-v_e)
+sC_{Q,c,r}v_c+sC_{Q,b,r}v_b+sC_{Q,e,r}v_e
$$

Compact matrix form:

$$
\boldsymbol{i}_{Q}=G_Q\boldsymbol{u}_{Q}+sC_Q\boldsymbol{v}_{Q}
$$

where the vectors are defined as follows:

$$
\boldsymbol{i}_{Q}=\left(i_{Q,c},i_{Q,b},i_{Q,e}\right)^{T}
$$

$$
\boldsymbol{u}_{Q}=\left(v_b-v_e,\ v_c-v_e\right)^{T}
$$

$$
\boldsymbol{v}_{Q}=\left(v_c,v_b,v_e\right)^{T}
$$

The rows of $G_Q$ are ordered as $c,b,e$, and its columns are ordered as $be,ce$:

| row | $be$ | $ce$ |
|---|---|---|
| $c$ | $g_{Q,be,c}$ | $g_{Q,ce,c}$ |
| $b$ | $g_{Q,be,b}$ | $g_{Q,ce,b}$ |
| $e$ | $g_{Q,be,e}$ | $g_{Q,ce,e}$ |

The rows of $C_Q$ are ordered as $c,b,e$, and its columns are ordered as $c,b,e$:

| row | $c$ | $b$ | $e$ |
|---|---|---|---|
| $c$ | $C_{Q,c,c}$ | $C_{Q,b,c}$ | $C_{Q,e,c}$ |
| $b$ | $C_{Q,c,b}$ | $C_{Q,b,b}$ | $C_{Q,e,b}$ |
| $e$ | $C_{Q,c,e}$ | $C_{Q,b,e}$ | $C_{Q,e,e}$ |

### 7.2 Equivalent BJT: `eq_circuit=True`

The equivalent-circuit stamp contains:

- a controlled current source from collector to emitter,
- an output resistance between collector and emitter,
- an input resistance between base and emitter,
- capacitances between base-emitter, base-collector, and collector-emitter.

Conductance/current terms:

$$
i_{gm}=g_Q(v_b-v_e),\qquad i_{gm}:c\rightarrow e
$$

$$
i_{ro}=\frac{v_c-v_e}{r_{Q,o}},\qquad i_{ro}:c\rightarrow e
$$

$$
i_{\pi}=\frac{v_b-v_e}{r_{Q,\pi}},\qquad i_{\pi}:b\rightarrow e
$$

Capacitance pairs:

$$
(b,e),\quad (b,c),\quad (c,e)
$$

Terminal current contributions:

$$
i_c=g_Q(v_b-v_e)+\frac{v_c-v_e}{r_{Q,o}}
+sC_{Q,bc}(v_c-v_b)+sC_{Q,ce}(v_c-v_e)
$$

$$
i_b=\frac{v_b-v_e}{r_{Q,\pi}}
+sC_{Q,be}(v_b-v_e)+sC_{Q,bc}(v_b-v_c)
$$

$$
i_e=-g_Q(v_b-v_e)-\frac{v_c-v_e}{r_{Q,o}}+\frac{v_e-v_b}{r_{Q,\pi}}
+sC_{Q,be}(v_e-v_b)+sC_{Q,ce}(v_e-v_c)
$$

## 8. Device-Off Options

`SSAS.solve(...)` accepts `off_device` and `off_device_cap`. These options apply to compact-device instances. They are not general netlist deletion controls.

Device names are normalized through the same alias rules described in `USAGE.md`. For example, `mn1` and `M_n1` identify the same MOSFET.

The two options remove different groups of the stamp:

| Options | Conductance/current group | Capacitance group |
|---|---|---|
| neither option | included | included |
| `off_device` contains device | removed | included |
| `off_device_cap` contains device | included | removed |
| both contain device | removed | removed |

For a single-device common-source stage, using both options on that MOSFET removes every MOSFET contribution. The remaining response can therefore reduce to zero even though the netlist instance is still present syntactically.

### 8.1 `off_device`

For selected MOSFETs in coefficient mode, `off_device` removes the conductance terms:

$$
g_{M,gs,r}(v_g-v_s),\quad
g_{M,ds,r}(v_d-v_s),\quad
g_{M,bs,r}(v_b-v_s)
$$

It leaves the MOSFET coefficient-mode device differential capacitance terms:

$$
sC_{M,c,r}v_c
$$

For selected MOSFETs in equivalent-circuit mode, `off_device` removes:

$$
g_M(v_g-v_s),\quad g_{M,mb}(v_b-v_s),\quad \frac{v_d-v_s}{r_M}
$$

It leaves all enabled MOSFET pair capacitances:

$$
C_{M,gs},\quad C_{M,gd},\quad C_{M,ds},\quad C_{M,gb},\quad C_{M,db},\quad C_{M,sb}
$$

For selected BJTs in `eq_circuit=False` mode, `off_device` removes:

$$
g_{Q,be,r}(v_b-v_e),\quad g_{Q,ce,r}(v_c-v_e)
$$

It leaves the BJT device differential capacitance terms:

$$
sC_{Q,k,r}v_k
$$

For selected BJTs in equivalent-circuit mode, `off_device` removes:

$$
g_Q(v_b-v_e),\quad \frac{v_c-v_e}{r_{Q,o}},\quad \frac{v_b-v_e}{r_{Q,\pi}}
$$

It leaves:

$$
C_{Q,be},\quad C_{Q,bc},\quad C_{Q,ce}
$$


### 8.2 `off_device_cap`

For selected MOSFETs in coefficient mode, `off_device_cap` removes:

$$
sC_{M,c,r}v_c
$$

It leaves:

$$
g_{M,gs,r}(v_g-v_s),\quad
g_{M,ds,r}(v_d-v_s),\quad
g_{M,bs,r}(v_b-v_s)
$$

For selected MOSFETs in equivalent-circuit mode, `off_device_cap` removes all enabled pair capacitances:

$$
C_{M,gs},\quad C_{M,gd},\quad C_{M,ds},\quad C_{M,gb},\quad C_{M,db},\quad C_{M,sb}
$$

It leaves:

$$
g_M(v_g-v_s),\quad g_{M,mb}(v_b-v_s),\quad \frac{v_d-v_s}{r_M}
$$

For selected BJTs in `eq_circuit=False` mode, `off_device_cap` removes:

$$
sC_{Q,k,r}v_k
$$

It leaves:

$$
g_{Q,be,r}(v_b-v_e),\quad g_{Q,ce,r}(v_c-v_e)
$$

For selected BJTs in equivalent-circuit mode, `off_device_cap` removes:

$$
C_{Q,be},\quad C_{Q,bc},\quad C_{Q,ce}
$$

It leaves:

$$
g_Q(v_b-v_e),\quad \frac{v_c-v_e}{r_{Q,o}},\quad \frac{v_b-v_e}{r_{Q,\pi}}
$$


## 9. Zero-Element Options

`zero_cap`, `zero_resistance`, and `zero_inductance` set selected symbolic element values to zero with their circuit meaning preserved.

For a selected passive or equivalent-circuit capacitance, zero means open circuit. For a coefficient-mode device capacitance term, zero means coefficient isolation, not a topology change. In both cases, the selected current contribution is not stamped:

$$
C_{M,gd}\in\mathrm{zero\_cap}\quad\Rightarrow\quad
sC_{M,gd}(v_g-v_d)\ \mathrm{is\ omitted}
$$


For a selected resistance, zero means ideal short circuit. The reciprocal admittance stamp is replaced by a branch-current unknown and a voltage equality constraint:

$$
R_D\in\mathrm{zero\_resistance}\quad\Rightarrow\quad
\begin{cases}
i_{R_D}\ \mathrm{is\ stamped\ into\ KCL}\\
v_a-v_b=0
\end{cases}
$$

$$
r_M\in\mathrm{zero\_resistance}\quad\Rightarrow\quad
\begin{cases}
i_{M,r}\ \mathrm{is\ stamped\ into\ KCL}\\
v_d-v_s=0
\end{cases}
$$

For a selected inductance, zero also means ideal short circuit. The branch current is still stamped, but the `sL` term is removed from the inductor constraint:

$$
L_D\in\mathrm{zero\_inductance}\quad\Rightarrow\quad
\begin{cases}
i_{L_D}\ \mathrm{is\ stamped\ into\ KCL}\\
v_a-v_b=0
\end{cases}
$$

Thus `zero_cap` changes topology for passive capacitors and equivalent-circuit terminal capacitors. In coefficient mode it removes selected device differential capacitance terms. `zero_resistance` and `zero_inductance` change topology by shorting selected passive elements or equivalent-circuit resistances.

## 10. Equal-Passive Option

`SSAS.solve(...)` also accepts `equal_passive`. This option maps multiple passive resistor, capacitor, or inductor value symbols to one alias only in the final rendered expression. The MNA stamps, topology, zero-element handling, and public current probes remain unchanged.

For example, `equal_passive=("R_X", ["r1", "r2"])` causes the final expression to substitute the resistor value symbols of `r1` and `r2` with `R_X` before TeX simplification. The option accepts one `(alias, equal_device)` pair, an iterable of such pairs, or a dictionary group with `{"alias": ..., "equal_device": ...}`.

Each alias group must contain only one passive type. Resistors, capacitors, and inductors cannot be mixed in one group, and non-passive devices are rejected. Zero-element options still refer to the original passive symbols because `equal_passive` is a final-output substitution.

## 11. Current Probe Support

An `SSAS` response request using `I(device)` is supported only for:

$$
V,\quad I,\quad R,\quad C,\quad L,\quad D
$$

It is not public for:

$$
M,\quad Q,\quad E,\quad G,\quad F,\quad H
$$

MOSFET, bipolar junction transistor, and controlled-source currents can affect the MNA system, but they are not exposed as public `I(device)` probes. Diode current is exposed because the diode primitive has a two-terminal current convention.

Ideal switch primitives are not supported. SPICE switch-style `S` and `W` devices do not have stamps in this extractor.
