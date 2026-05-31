from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
import math
import re
from typing import Iterable


ZERO = 0
ONE = 1
RATIONAL = 2
ATOM = 3
NEG = 4
ADD = 5
MUL = 6
POW_INT = 7


class KernelError(Exception):
    status = "kernel_error"

    def __init__(self, message: str, status: str | None = None):
        super().__init__(message)
        self.message = message
        if status is not None:
            self.status = status


class KernelExactDivisionError(KernelError):
    status = "exact_division_failed"


class KernelSolveError(KernelError):
    status = "singular_or_no_closed_form"


class KernelRequestError(KernelError):
    status = "invalid_request"


class KernelTexError(KernelError):
    status = "tex_export_failed"


class KernelZeroDenominatorError(KernelError):
    status = "zero_denominator"


def plain(text: str) -> str:
    return re.sub(r"[^0-9A-Za-z]+", "_", str(text)).strip("_")


def nodeName(node: str) -> str:
    return plain(node).lower()


def splitTopLevelSlash(text: str) -> tuple[str, str] | None:
    depth = 0
    for idx, char in enumerate(text):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif char == "/" and depth == 0:
            return text[:idx].strip(), text[idx + 1 :].strip()
    return None


@dataclass
class AtomTable:
    name_to_id: dict[str, int] = field(default_factory=dict)
    names: list[str] = field(default_factory=list)
    tex_names: list[str] = field(default_factory=list)
    kind: list[str] = field(default_factory=list)
    base_param: list[int | None] = field(default_factory=list)

    def param(self, name: str) -> int:
        key = plain(name) or str(name)
        existing = self.name_to_id.get(key)
        if existing is not None:
            return existing
        idx = len(self.names)
        self.name_to_id[key] = idx
        self.names.append(key)
        self.tex_names.append(texName(name))
        self.kind.append("param")
        self.base_param.append(None)
        return idx

    def invParam(self, name: str) -> int:
        base = self.param(name)
        key = f"inv_{self.names[base]}"
        existing = self.name_to_id.get(key)
        if existing is not None:
            return existing
        idx = len(self.names)
        self.name_to_id[key] = idx
        self.names.append(key)
        self.tex_names.append(r"\frac{1}{" + self.tex_names[base] + "}")
        self.kind.append("reciprocal_atom")
        self.base_param.append(base)
        return idx

    def atomLatex(self, atom_id: int) -> str:
        return self.tex_names[atom_id]


def texName(name: str) -> str:
    text = str(name).strip()
    if "{" in text or "\\" in text:
        return text
    key = plain(text)
    if "_" not in key:
        return key
    head, tail = key.split("_", 1)
    return f"{head}_{{{tail}}}"


class ExprArena:
    # Expressions are interned by structural key, so an expr_id is also a dense
    # index into immutable node metadata.  ADD/MUL nodes are always built by
    # their canonical constructors; callers can rely on sorted children and a
    # leading rational coefficient for products when one exists.
    def __init__(self, atoms: AtomTable):
        self.atoms = atoms
        self.kind: list[int] = []
        self.value: list[object] = []
        self.child_tuple: list[tuple[int, ...]] = []
        self._cache: dict[tuple, int] = {}
        self._sort_key_cache: list[tuple | None] = []
        self._factor_power_items_cache: list[tuple[tuple[int, int], ...] | None] = []
        self._common_term_factor_items_cache: list[tuple[tuple[int, int], ...] | None] = []
        self._common_atom_cache: list[frozenset[int] | None] = []
        self._binary_add_cache: dict[tuple[int, int], int] = {}
        self._binary_mul_cache: dict[tuple[int, int], int] = {}
        self.zero = self.internKey((ZERO, 0, ()))
        self.one = self.internKey((ONE, 1, ()))

    def internKey(self, key: tuple) -> int:
        existing = self._cache.get(key)
        if existing is not None:
            return existing
        idx = len(self.kind)
        self._cache[key] = idx
        self.kind.append(int(key[0]))
        self.value.append(key[1])
        child_tuple = tuple(key[2])
        self.child_tuple.append(child_tuple)
        self._sort_key_cache.append(None)
        self._factor_power_items_cache.append(None)
        self._common_term_factor_items_cache.append(None)
        self._common_atom_cache.append(None)
        return idx

    def nodeChildren(self, expr_id: int) -> tuple[int, ...]:
        return self.child_tuple[expr_id]

    def rational(self, value: int | Fraction | float | str) -> int:
        frac = parseFraction(value)
        if frac == 0:
            return self.zero
        if frac == 1:
            return self.one
        return self.internKey((RATIONAL, (frac.numerator, frac.denominator), ()))

    def atom(self, atom_id: int) -> int:
        return self.internKey((ATOM, int(atom_id), ()))

    def neg(self, expr_id: int) -> int:
        return self.mul(self.rational(-1), expr_id)

    def sub(self, left: int, right: int) -> int:
        return self.add(left, self.neg(right))

    def add(self, *expr_ids: int) -> int:
        if len(expr_ids) == 2:
            left, right = expr_ids
            if left == self.zero:
                return right
            if right == self.zero:
                return left
            key = (left, right) if left <= right else (right, left)
            existing = self._binary_add_cache.get(key)
            if existing is not None:
                return existing
            result = self._add_uncached(left, right)
            self._binary_add_cache[key] = result
            return result
        return self._add_uncached(*expr_ids)

    def _add_uncached(self, *expr_ids: int) -> int:
        terms: list[int] = []
        kind = self.kind
        node_children = self.nodeChildren
        for expr_id in expr_ids:
            if expr_id == self.zero:
                continue
            if kind[expr_id] == ADD:
                terms.extend(node_children(expr_id))
            else:
                terms.append(expr_id)
        if not terms:
            return self.zero
        collected: dict[int, Fraction] = {}
        for term in terms:
            coeff, base = self.splitCoeff(term)
            collected[base] = collected.get(base, Fraction(0)) + coeff
        out: list[int] = []
        for base, coeff in collected.items():
            if coeff == 0:
                continue
            if base == self.one:
                out.append(self.rational(coeff))
            elif coeff == 1:
                out.append(base)
            else:
                out.append(self.mul(self.rational(coeff), base))
        if not out:
            return self.zero
        out.sort(key=self.sortKey)
        if len(out) == 1:
            return out[0]
        return self.internKey((ADD, None, tuple(out)))

    def mul(self, *expr_ids: int) -> int:
        if len(expr_ids) == 2:
            left, right = expr_ids
            if left == self.zero or right == self.zero:
                return self.zero
            if left == self.one:
                return right
            if right == self.one:
                return left
            key = (left, right) if left <= right else (right, left)
            existing = self._binary_mul_cache.get(key)
            if existing is not None:
                return existing
            result = self._mul_uncached(left, right)
            self._binary_mul_cache[key] = result
            return result
        return self._mul_uncached(*expr_ids)

    def _mul_uncached(self, *expr_ids: int) -> int:
        factors: list[int] = []
        rational = Fraction(1)
        kind = self.kind
        value = self.value
        node_children = self.nodeChildren
        for expr_id in expr_ids:
            if expr_id == self.zero:
                return self.zero
            if expr_id == self.one:
                continue
            expr_kind = kind[expr_id]
            if expr_kind == RATIONAL:
                p, q = value[expr_id]
                rational *= Fraction(int(p), int(q))
            elif expr_kind == MUL:
                factors.extend(node_children(expr_id))
            else:
                factors.append(expr_id)
        if rational == 0:
            return self.zero
        power_map: dict[int, int] = {}
        for factor in factors:
            if factor == self.one:
                continue
            factor_kind = kind[factor]
            if factor_kind == RATIONAL:
                p, q = value[factor]
                rational *= Fraction(int(p), int(q))
                continue
            if factor_kind == POW_INT:
                base, exp = value[factor]
                power_map[int(base)] = power_map.get(int(base), 0) + int(exp)
            else:
                power_map[factor] = power_map.get(factor, 0) + 1
        out: list[int] = []
        if rational != 1:
            out.append(self.rational(rational))
        for base, exp in power_map.items():
            if exp == 0:
                continue
            out.append(base if exp == 1 else self.powInt(base, exp))
        if not out:
            return self.one
        out.sort(key=self.sortKey)
        if len(out) == 1:
            return out[0]
        return self.internKey((MUL, None, tuple(out)))

    def powInt(self, base: int, exponent: int) -> int:
        if exponent < 0:
            raise KernelError("POW_INT accepts only nonnegative exponents")
        if exponent == 0:
            return self.one
        if exponent == 1:
            return base
        if base == self.zero:
            return self.zero
        if base == self.one:
            return self.one
        return self.internKey((POW_INT, (base, int(exponent)), (base,)))

    def splitCoeff(self, expr_id: int) -> tuple[Fraction, int]:
        if expr_id == self.one:
            return Fraction(1), self.one
        kind = self.kind
        value = self.value
        expr_kind = kind[expr_id]
        if expr_kind == RATIONAL:
            p, q = value[expr_id]
            return Fraction(int(p), int(q)), self.one
        if expr_kind == MUL:
            children = self.child_tuple[expr_id]
            first = children[0]
            if kind[first] == RATIONAL:
                p, q = value[first]
                rest = children[1:]
                if not rest:
                    return Fraction(int(p), int(q)), self.one
                if len(rest) == 1:
                    return Fraction(int(p), int(q)), rest[0]
                return Fraction(int(p), int(q)), self.internKey((MUL, None, rest))
        return Fraction(1), expr_id

    def sortKey(self, expr_id: int) -> tuple:
        # Sort keys are part of the rendered text contract.  The cached value is
        # exactly the recursive tuple that the uncached implementation returned.
        cached = self._sort_key_cache[expr_id]
        if cached is not None:
            return cached
        kind = self.kind[expr_id]
        value = self.value[expr_id]
        if kind in {ZERO, ONE, RATIONAL, ATOM}:
            out = (kind, value)
        elif kind == POW_INT:
            base, exp = value
            out = (kind, exp, self.sortKey(base))
        else:
            out = (kind, tuple(self.sortKey(child) for child in self.child_tuple[expr_id]))
        self._sort_key_cache[expr_id] = out
        return out

    def commonAtomExprs(self, expr_id: int) -> frozenset[int]:
        cached = self._common_atom_cache[expr_id]
        if cached is not None:
            return cached
        if expr_id == self.zero:
            out = frozenset()
            self._common_atom_cache[expr_id] = out
            return out
        kind = self.kind[expr_id]
        if kind == ATOM:
            out = frozenset((expr_id,))
        elif kind == POW_INT:
            base, exp = self.value[expr_id]
            base = int(base)
            out = frozenset((base,)) if exp > 0 and self.kind[base] == ATOM else self.commonAtomExprs(base)
        elif kind == MUL:
            out: set[int] = set()
            for child in self.child_tuple[expr_id]:
                out.update(self.commonAtomExprs(child))
            out = frozenset(out)
        elif kind == ADD:
            children = self.child_tuple[expr_id]
            common: frozenset[int] | None = None
            for child in children:
                if child == self.zero:
                    continue
                child_atoms = self.commonAtomExprs(child)
                common = child_atoms if common is None else common & child_atoms
                if not common:
                    break
            out = common if common is not None else frozenset()
        else:
            out = frozenset()
        self._common_atom_cache[expr_id] = out
        return out

    def divByAtomExpr(self, expr_id: int, atom_expr: int) -> int | None:
        if expr_id == self.zero:
            return self.zero
        if expr_id == atom_expr:
            return self.one
        kind = self.kind[expr_id]
        if kind == POW_INT:
            base, exp = self.value[expr_id]
            if int(base) == atom_expr and int(exp) > 0:
                return self.powInt(atom_expr, int(exp) - 1)
        if kind == MUL:
            children = list(self.nodeChildren(expr_id))
            for idx, child in enumerate(children):
                if child == atom_expr:
                    return self.mul(*(children[:idx] + children[idx + 1 :]))
                if self.kind[child] == POW_INT:
                    base, exp = self.value[child]
                    if int(base) == atom_expr and int(exp) > 0:
                        replacement = self.powInt(atom_expr, int(exp) - 1)
                        children[idx] = replacement
                        return self.mul(*children)
            for idx, child in enumerate(children):
                quotient = self.divByAtomExpr(child, atom_expr)
                if quotient is not None:
                    children[idx] = quotient
                    return self.mul(*children)
        if kind == ADD:
            quotients: list[int] = []
            for child in self.nodeChildren(expr_id):
                quotient = self.divByAtomExpr(child, atom_expr)
                if quotient is None:
                    return None
                quotients.append(quotient)
            return self.add(*quotients)
        return None

    def divByExpr(self, expr_id: int, divisor: int) -> int | None:
        if divisor == self.one:
            return expr_id
        if expr_id == self.zero:
            return self.zero
        if expr_id == divisor:
            return self.one
        divisor_kind = self.kind[divisor]
        if divisor_kind == RATIONAL:
            p, q = self.value[divisor]
            if p == 0:
                return None
            if self.kind[expr_id] == ADD:
                quotients: list[int] = []
                for child in self.nodeChildren(expr_id):
                    quotient = self.divByExpr(child, divisor)
                    if quotient is None:
                        return None
                    quotients.append(quotient)
                return self.add(*quotients)
            return self.mul(expr_id, self.rational(Fraction(int(q), int(p))))
        if divisor_kind == ATOM:
            return self.divByAtomExpr(expr_id, divisor)
        if divisor_kind == MUL:
            result = expr_id
            for factor in self.nodeChildren(divisor):
                result = self.divByExpr(result, factor)
                if result is None:
                    return None
            return result
        if divisor_kind == POW_INT:
            base, exp = self.value[divisor]
            result = expr_id
            for _ in range(int(exp)):
                result = self.divByExpr(result, int(base))
                if result is None:
                    return None
            return result
        kind = self.kind[expr_id]
        if kind == MUL:
            children = list(self.nodeChildren(expr_id))
            for idx, child in enumerate(children):
                if child == divisor:
                    return self.mul(*(children[:idx] + children[idx + 1 :]))
                if self.kind[child] == POW_INT and self.kind[divisor] == POW_INT:
                    base, exp = self.value[child]
                    div_base, div_exp = self.value[divisor]
                    if int(base) == int(div_base) and int(exp) >= int(div_exp):
                        replacement = self.powInt(int(base), int(exp) - int(div_exp))
                        children[idx] = replacement
                        return self.mul(*children)
                quotient = self.divByExpr(child, divisor)
                if quotient is not None:
                    children[idx] = quotient
                    return self.mul(*children)
        if kind == ADD:
            quotients: list[int] = []
            for child in self.nodeChildren(expr_id):
                quotient = self.divByExpr(child, divisor)
                if quotient is None:
                    return None
                quotients.append(quotient)
            return self.add(*quotients)
        return None

    def factorPowerItems(self, expr_id: int) -> tuple[tuple[int, int], ...]:
        cached = self._factor_power_items_cache[expr_id]
        if cached is not None:
            return cached
        if expr_id in (self.zero, self.one):
            out: tuple[tuple[int, int], ...] = ()
            self._factor_power_items_cache[expr_id] = out
            return out
        kind = self.kind[expr_id]
        if kind == RATIONAL:
            out = ()
            self._factor_power_items_cache[expr_id] = out
            return out
        if kind == MUL:
            out: dict[int, int] = {}
            for child in self.child_tuple[expr_id]:
                for factor, power in self.factorPowerItems(child):
                    out[factor] = out.get(factor, 0) + power
            items = tuple(out.items())
        elif kind == POW_INT:
            base, exp = self.value[expr_id]
            if int(exp) <= 0:
                items = ()
            else:
                items = ((int(base), int(exp)),)
        else:
            items = ((expr_id, 1),)
        self._factor_power_items_cache[expr_id] = items
        return items

    def commonTermFactorPowerItems(self, expr_id: int) -> tuple[tuple[int, int], ...]:
        # Non-ADD expressions use the normal factor-power view.  ADD nodes keep
        # only factors present in every nonzero term.
        cached = self._common_term_factor_items_cache[expr_id]
        if cached is not None:
            return cached
        if self.kind[expr_id] != ADD:
            out = self.factorPowerItems(expr_id)
            self._common_term_factor_items_cache[expr_id] = out
            return out
        children = self.child_tuple[expr_id]
        common: dict[int, int] | None = None
        for child in children:
            if child == self.zero:
                continue
            current_items = self.factorPowerItems(child)
            if common is None:
                common = dict(current_items)
                continue
            current = dict(current_items)
            for factor in list(common):
                if factor not in current:
                    del common[factor]
                else:
                    common[factor] = min(common[factor], current[factor])
            if not common:
                break
        out = tuple(common.items()) if common is not None else ()
        self._common_term_factor_items_cache[expr_id] = out
        return out

    def isNegative(self, expr_id: int) -> bool:
        kind = self.kind
        value = self.value
        if kind[expr_id] == RATIONAL:
            p, _ = value[expr_id]
            return int(p) < 0
        if kind[expr_id] == MUL:
            first = self.child_tuple[expr_id][0]
            if kind[first] == RATIONAL:
                p, _ = value[first]
                return int(p) < 0
        return False

    def evaluate(self, expr_id: int, values: dict[str, Fraction]) -> Fraction:
        kind = self.kind[expr_id]
        if kind == ZERO:
            return Fraction(0)
        if kind == ONE:
            return Fraction(1)
        if kind == RATIONAL:
            p, q = self.value[expr_id]
            return Fraction(int(p), int(q))
        if kind == ATOM:
            atom_id = int(self.value[expr_id])
            if self.atoms.kind[atom_id] == "reciprocal_atom":
                base = self.atoms.base_param[atom_id]
                if base is None:
                    raise KernelError("reciprocal atom has no base parameter")
                base_name = self.atoms.names[base]
                if base_name not in values:
                    raise KernelError(f"missing numeric value for {base_name}")
                return Fraction(1, 1) / values[base_name]
            name = self.atoms.names[atom_id]
            if name not in values:
                raise KernelError(f"missing numeric value for {name}")
            return values[name]
        if kind == ADD:
            total = Fraction(0)
            for child in self.nodeChildren(expr_id):
                total += self.evaluate(child, values)
            return total
        if kind == MUL:
            result = Fraction(1)
            for child in self.nodeChildren(expr_id):
                result *= self.evaluate(child, values)
            return result
        if kind == POW_INT:
            base, exp = self.value[expr_id]
            return self.evaluate(int(base), values) ** int(exp)
        raise KernelError("unknown expression node")


def parseFraction(value: int | Fraction | float | str) -> Fraction:
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value, 1)
    if isinstance(value, float):
        return Fraction(str(value))
    text = str(value).strip()
    if text.startswith("+"):
        text = text[1:]
    return Fraction(text)


class SPolyArena:
    # Polynomials are also dense-id interned.  Operation caches stay keyed by
    # polynomial ids; structural analysis caches use id-indexed lists.
    def __init__(self, expr: ExprArena):
        self.expr = expr
        self.coeffs: list[tuple[int, ...]] = []
        self._cache: dict[tuple[int, ...], int] = {}
        self._add_cache: dict[tuple[int, int], int] = {}
        self._mul_cache: dict[tuple[int, int], int] = {}
        self._neg_cache: list[int | None] = []
        self._common_atom_cache: list[frozenset[int] | None] = []
        self.known_factor_children: dict[int, tuple[int, ...]] = {}
        self.zero = self.poly(())
        self.one = self.poly((expr.one,))

    def poly(self, coeffs: Iterable[int]) -> int:
        data = list(coeffs)
        while data and data[-1] == self.expr.zero:
            data.pop()
        key = tuple(data)
        existing = self._cache.get(key)
        if existing is not None:
            return existing
        idx = len(self.coeffs)
        self._cache[key] = idx
        self.coeffs.append(key)
        self._neg_cache.append(None)
        self._common_atom_cache.append(None)
        return idx

    def fromExpr(self, expr_id: int) -> int:
        return self.zero if expr_id == self.expr.zero else self.poly((expr_id,))

    def degree(self, poly_id: int) -> int:
        return len(self.coeffs[poly_id]) - 1

    def add(self, left: int, right: int) -> int:
        if left == self.zero:
            return right
        if right == self.zero:
            return left
        key = (left, right) if left <= right else (right, left)
        existing = self._add_cache.get(key)
        if existing is not None:
            return existing
        a = self.coeffs[left]
        b = self.coeffs[right]
        n = max(len(a), len(b))
        out = []
        for idx in range(n):
            ca = a[idx] if idx < len(a) else self.expr.zero
            cb = b[idx] if idx < len(b) else self.expr.zero
            out.append(self.expr.add(ca, cb))
        result = self.poly(out)
        self._add_cache[key] = result
        return result

    def neg(self, poly_id: int) -> int:
        if poly_id == self.zero:
            return self.zero
        cached = self._neg_cache[poly_id]
        if cached is not None:
            return cached
        result = self.poly(self.expr.neg(coeff) for coeff in self.coeffs[poly_id])
        self._neg_cache[poly_id] = result
        return result

    def sub(self, left: int, right: int) -> int:
        return self.add(left, self.neg(right))

    def mul(self, left: int, right: int) -> int:
        if left == self.zero or right == self.zero:
            return self.zero
        if left == self.one:
            return right
        if right == self.one:
            return left
        key = (left, right) if left <= right else (right, left)
        existing = self._mul_cache.get(key)
        if existing is not None:
            return existing
        a = self.coeffs[left]
        b = self.coeffs[right]
        out = [self.expr.zero] * (len(a) + len(b) - 1)
        for i, ca in enumerate(a):
            if ca == self.expr.zero:
                continue
            for j, cb in enumerate(b):
                if cb == self.expr.zero:
                    continue
                out[i + j] = self.expr.add(out[i + j], self.expr.mul(ca, cb))
        result = self.poly(out)
        self.rememberProductFactors(result, left, right)
        self._mul_cache[key] = result
        return result

    def scale(self, expr_id: int, poly_id: int) -> int:
        if expr_id == self.expr.zero or poly_id == self.zero:
            return self.zero
        if expr_id == self.expr.one:
            return poly_id
        result = self.poly(self.expr.mul(expr_id, coeff) for coeff in self.coeffs[poly_id])
        if self.degree(poly_id) > 0:
            self.known_factor_children[result] = self.known_factor_children.get(poly_id, (poly_id,))
        return result

    def shift(self, poly_id: int, amount: int) -> int:
        if poly_id == self.zero:
            return self.zero
        result = self.poly([self.expr.zero] * amount + list(self.coeffs[poly_id]))
        if amount > 0:
            s_factor = self.poly((self.expr.zero, self.expr.one))
            factors = tuple([s_factor] * amount + list(self.known_factor_children.get(poly_id, (poly_id,))))
            self.known_factor_children[result] = factors
        return result

    def leadingCoeff(self, poly_id: int) -> int:
        if poly_id == self.zero:
            return self.expr.zero
        return self.coeffs[poly_id][-1]

    def exactDiv(self, numerator: int, denominator: int) -> tuple[int, bool]:
        if denominator == self.zero:
            raise KernelExactDivisionError("division by zero polynomial")
        if numerator == self.zero:
            return self.zero, True
        if denominator == self.one:
            return numerator, True
        if numerator == denominator:
            return self.one, True
        den_coeffs = self.coeffs[denominator]
        if len(den_coeffs) == 1:
            divisor = den_coeffs[0]
            out: list[int] = []
            for coeff in self.coeffs[numerator]:
                quotient = self.expr.divByExpr(coeff, divisor)
                if quotient is None:
                    return self.zero, False
                out.append(quotient)
            return self.poly(out), True
        rem = list(self.coeffs[numerator])
        den_degree = len(den_coeffs) - 1
        den_lead = den_coeffs[-1]
        if len(rem) < len(den_coeffs):
            return self.zero, False
        quotient = [self.expr.zero] * (len(rem) - den_degree)
        while len(rem) >= len(den_coeffs):
            rem_degree = len(rem) - 1
            lead = rem[-1]
            if lead == self.expr.zero:
                rem.pop()
                continue
            q_coeff = self.expr.divByExpr(lead, den_lead)
            if q_coeff is None:
                return self.zero, False
            shift = rem_degree - den_degree
            quotient[shift] = self.expr.add(quotient[shift], q_coeff)
            for idx, den_coeff in enumerate(den_coeffs):
                pos = idx + shift
                rem[pos] = self.expr.sub(rem[pos], self.expr.mul(q_coeff, den_coeff))
            while rem and rem[-1] == self.expr.zero:
                rem.pop()
        if any(coeff != self.expr.zero for coeff in rem):
            return self.zero, False
        return self.poly(quotient), True

    def commonAtomExprs(self, poly_id: int) -> frozenset[int]:
        cached = self._common_atom_cache[poly_id]
        if cached is not None:
            return cached
        common: frozenset[int] | None = None
        for coeff in self.coeffs[poly_id]:
            if coeff == self.expr.zero:
                continue
            coeff_atoms = self.expr.commonAtomExprs(coeff)
            common = coeff_atoms if common is None else common & coeff_atoms
            if not common:
                break
        out = common if common is not None else frozenset()
        self._common_atom_cache[poly_id] = out
        return out

    def divByAtomExpr(self, poly_id: int, atom_expr: int) -> int | None:
        out: list[int] = []
        for coeff in self.coeffs[poly_id]:
            quotient = self.expr.divByAtomExpr(coeff, atom_expr)
            if quotient is None:
                return None
            out.append(quotient)
        return self.poly(out)

    def factorChildrenFor(self, poly_id: int) -> tuple[int, ...]:
        children = self.known_factor_children.get(poly_id)
        if children is None:
            return (poly_id,)
        return children

    def rememberProductFactors(self, result: int, left: int, right: int) -> None:
        if result == self.zero:
            return
        factors: list[int] = []
        for poly_id in (left, right):
            if self.degree(poly_id) <= 0:
                continue
            factors.extend(self.factorChildrenFor(poly_id))
        if len(factors) >= 2:
            self.known_factor_children[result] = tuple(factors)

    def evaluate(self, poly_id: int, values: dict[str, Fraction], s_value: Fraction) -> Fraction:
        result = Fraction(0)
        power = Fraction(1)
        for coeff in self.coeffs[poly_id]:
            result += self.expr.evaluate(coeff, values) * power
            power *= s_value
        return result


@dataclass(frozen=True)
class SRat:
    num: int
    den: int


@dataclass
class UnknownInfo:
    name: str
    kind: str
    node: str = ""
    device: str = ""


@dataclass
class Response:
    request: str
    value: SRat
    kernel: "SolverKernel"
    passive_value_aliases: dict[str, str] = field(default_factory=dict)

    def evaluate(self, assignments: dict[str, int | Fraction], s_value: int | Fraction) -> Fraction:
        values = {plain(key): parseFraction(value) for key, value in assignments.items()}
        return self.kernel.evaluateRat(self.value, values, parseFraction(s_value))


class SolverKernel:
    def __init__(self, circuit):
        self.circuit = circuit
        self.atoms: AtomTable = circuit.atoms
        self.expr: ExprArena = circuit.expr
        self.spoly: SPolyArena = circuit.spoly
        self.unknowns: list[UnknownInfo] = circuit.unknowns
        self.equations = circuit.equations
        self.solution: list[SRat] | None = None
        self._zero_rat = SRat(self.spoly.zero, self.spoly.one)
        self._one_rat = SRat(self.spoly.one, self.spoly.one)
        self._det_cache: dict[tuple[tuple[int, ...], ...], SRat] = {}
        self._reduced_system = None

    def ratFromPoly(self, poly_id: int) -> SRat:
        return SRat(poly_id, self.spoly.one)

    def ratNeg(self, value: SRat) -> SRat:
        return self.normalizeRat(SRat(self.spoly.neg(value.num), value.den))

    def ratAdd(self, left: SRat, right: SRat) -> SRat:
        if left.num == self.spoly.zero:
            return right
        if right.num == self.spoly.zero:
            return left
        if left.den == right.den:
            return self.normalizeRat(SRat(self.spoly.add(left.num, right.num), left.den))
        num = self.spoly.add(self.spoly.mul(left.num, right.den), self.spoly.mul(right.num, left.den))
        den = self.spoly.mul(left.den, right.den)
        return self.normalizeRat(SRat(num, den))

    def ratSub(self, left: SRat, right: SRat) -> SRat:
        return self.ratAdd(left, self.ratNeg(right))

    def ratMul(self, left: SRat, right: SRat) -> SRat:
        if left.num == self.spoly.zero or right.num == self.spoly.zero:
            return self._zero_rat
        return self.normalizeRat(SRat(self.spoly.mul(left.num, right.num), self.spoly.mul(left.den, right.den)))

    def ratDiv(self, left: SRat, right: SRat) -> SRat:
        if right.num == self.spoly.zero:
            raise KernelZeroDenominatorError("zero denominator in response request")
        if left.num == self.spoly.zero:
            return self._zero_rat
        if left.den == right.den:
            return self.normalizeRat(SRat(left.num, right.num))
        return self.normalizeRat(SRat(self.spoly.mul(left.num, right.den), self.spoly.mul(left.den, right.num)))

    def normalizeRat(self, value: SRat) -> SRat:
        if value.den == self.spoly.zero:
            raise KernelSolveError("zero denominator in rational expression")
        if value.num == self.spoly.zero:
            return self._zero_rat
        if value.num == value.den:
            return self._one_rat
        num = value.num
        den = value.den
        changed = True
        while changed:
            changed = False
            common = self.spoly.commonAtomExprs(num) & self.spoly.commonAtomExprs(den)
            for atom_expr in sorted(common):
                next_num = self.spoly.divByAtomExpr(num, atom_expr)
                next_den = self.spoly.divByAtomExpr(den, atom_expr)
                if next_num is not None and next_den is not None:
                    num, den = next_num, next_den
                    changed = True
                    break
        if num == den:
            return self._one_rat
        if den == self.spoly.one:
            return SRat(num, den)
        return SRat(num, den)

    def solve(self) -> list[SRat]:
        if self.solution is not None:
            return self.solution
        n = len(self.unknowns)
        equations, active, substitutions = self.substituteUnitEquations(self.equations, n)
        if len(equations) != len(active):
            raise KernelSolveError(f"reduced MNA system is not square: {len(equations)} equations for {len(active)} unknowns")
        active_list = sorted(active)
        active_pos = {unknown_id: idx for idx, unknown_id in enumerate(active_list)}
        matrix: list[list[SRat]] = []
        rhs: list[SRat] = []
        for equation in equations:
            row = [self._zero_rat for _ in active_list]
            for unknown_id, coeff in equation.coeffs.items():
                if unknown_id in active_pos:
                    row[active_pos[unknown_id]] = self.ratFromPoly(coeff)
                elif coeff != self.spoly.zero:
                    raise KernelSolveError("internal substitution left an inactive unknown in the reduced system")
            matrix.append(row)
            rhs.append(self.ratFromPoly(self.spoly.neg(equation.const)))
        reduced_solution = self.solveDense(matrix, rhs) if active_list else []
        solution = [self._zero_rat for _ in range(n)]
        for unknown_id, reduced_id in active_pos.items():
            solution[unknown_id] = reduced_solution[reduced_id]
        for unknown_id in sorted(substitutions):
            solution[unknown_id] = self.evaluateLinearWithSolution(substitutions[unknown_id], solution)
        self.solution = solution
        return solution

    def substituteUnitEquations(self, equations, unknown_count: int):
        work = [equation.copy() for equation in equations]
        active = set(range(unknown_count))
        substitutions: dict[int, object] = {}
        neg_one = self.spoly.neg(self.spoly.one)
        changed = True
        while changed:
            changed = False
            for row_id, equation in enumerate(work):
                candidates = [uid for uid, coeff in equation.coeffs.items() if uid in active and (coeff == self.spoly.one or coeff == neg_one)]
                if not candidates:
                    continue
                unknown_id = self.chooseSubstitutionUnknown(candidates)
                coeff = equation.coeffs[unknown_id]
                rest = equation.copy()
                rest.coeffs.pop(unknown_id, None)
                value = rest.neg() if coeff == self.spoly.one else rest
                work.pop(row_id)
                active.remove(unknown_id)
                for idx, other in enumerate(work):
                    work[idx] = self.substituteLinear(other, unknown_id, value)
                for key, old_value in list(substitutions.items()):
                    substitutions[key] = self.substituteLinear(old_value, unknown_id, value)
                substitutions[unknown_id] = value
                changed = True
                break
        return work, active, substitutions

    def chooseSubstitutionUnknown(self, candidates: list[int]) -> int:
        def score(unknown_id: int) -> tuple[int, int]:
            info = self.unknowns[unknown_id]
            branch_bonus = 0 if info.kind == "branch_current" else 1
            return (branch_bonus, unknown_id)
        return sorted(candidates, key=score)[0]

    def substituteLinear(self, equation, unknown_id: int, value):
        coeff = equation.coeffs.get(unknown_id)
        if coeff is None:
            return equation
        out = equation.copy()
        out.coeffs.pop(unknown_id, None)
        return out.add(value.scale(coeff))

    def evaluateLinearWithSolution(self, equation, solution: list[SRat]) -> SRat:
        total = self.ratFromPoly(equation.const)
        for unknown_id, coeff in equation.coeffs.items():
            total = self.ratAdd(total, self.ratMul(self.ratFromPoly(coeff), solution[unknown_id]))
        return total

    def solveDense(self, matrix: list[list[SRat]], rhs: list[SRat]) -> list[SRat]:
        n = len(matrix)
        col_perm = list(range(n))
        for k in range(n):
            pivot = self.choosePivot(matrix, k)
            if pivot is None:
                raise KernelSolveError("singular MNA system")
            pivot_row, pivot_col = pivot
            if pivot_row != k:
                matrix[k], matrix[pivot_row] = matrix[pivot_row], matrix[k]
                rhs[k], rhs[pivot_row] = rhs[pivot_row], rhs[k]
            if pivot_col != k:
                for row in matrix:
                    row[k], row[pivot_col] = row[pivot_col], row[k]
                col_perm[k], col_perm[pivot_col] = col_perm[pivot_col], col_perm[k]
            pivot_value = matrix[k][k]
            for row_id in range(k + 1, n):
                entry = matrix[row_id][k]
                if entry.num == self.spoly.zero:
                    continue
                factor = self.ratDiv(entry, pivot_value)
                matrix[row_id][k] = self._zero_rat
                for col_id in range(k + 1, n):
                    matrix[row_id][col_id] = self.ratSub(matrix[row_id][col_id], self.ratMul(factor, matrix[k][col_id]))
                rhs[row_id] = self.ratSub(rhs[row_id], self.ratMul(factor, rhs[k]))
        perm_solution = [self._zero_rat for _ in range(n)]
        for row_id in range(n - 1, -1, -1):
            total = rhs[row_id]
            for col_id in range(row_id + 1, n):
                if matrix[row_id][col_id].num != self.spoly.zero:
                    total = self.ratSub(total, self.ratMul(matrix[row_id][col_id], perm_solution[col_id]))
            perm_solution[row_id] = self.ratDiv(total, matrix[row_id][row_id])
        solution = [self._zero_rat for _ in range(n)]
        for col_id, original_id in enumerate(col_perm):
            solution[original_id] = perm_solution[col_id]
        return solution

    def choosePivot(self, matrix: list[list[SRat]], start: int) -> tuple[int, int] | None:
        best: tuple[tuple[int, int, int, int], int, int] | None = None
        n = len(matrix)
        for row_id in range(start, n):
            row_nnz = sum(1 for col_id in range(start, n) if matrix[row_id][col_id].num != self.spoly.zero)
            if row_nnz == 0:
                continue
            for col_id in range(start, n):
                value = matrix[row_id][col_id]
                if value.num == self.spoly.zero:
                    continue
                col_nnz = sum(1 for r in range(start, n) if matrix[r][col_id].num != self.spoly.zero)
                fill = (row_nnz - 1) * (col_nnz - 1)
                degree = max(0, self.spoly.degree(value.num)) + max(0, self.spoly.degree(value.den))
                size = len(self.spoly.coeffs[value.num]) + len(self.spoly.coeffs[value.den])
                score = (fill, degree, size, row_id + col_id)
                if best is None or score < best[0]:
                    best = (score, row_id, col_id)
        return None if best is None else (best[1], best[2])

    def responseFromProbes(self, request: str, numerator_probe, denominator_probe) -> Response:
        denominator = self.probeNumerator(denominator_probe)
        if denominator.num == self.spoly.zero:
            raise KernelZeroDenominatorError("zero denominator in response request")
        numerator = self.probeNumerator(numerator_probe)
        return Response(request, self.ratDiv(numerator, denominator), self)

    def prepareReducedSystem(self):
        if self._reduced_system is not None:
            return self._reduced_system
        equations, active, substitutions = self.substituteUnitEquations(self.equations, len(self.unknowns))
        if len(equations) != len(active):
            raise KernelSolveError(f"reduced MNA system is not square: {len(equations)} equations for {len(active)} unknowns")
        active_list = sorted(active)
        active_pos = {unknown_id: idx for idx, unknown_id in enumerate(active_list)}
        a_rows: list[list[int]] = []
        b_vec: list[int] = []
        for equation in equations:
            row = [self.spoly.zero for _ in active_list]
            for unknown_id, coeff in equation.coeffs.items():
                if unknown_id not in active_pos:
                    raise KernelSolveError("internal substitution left an inactive unknown in the reduced system")
                row[active_pos[unknown_id]] = coeff
            a_rows.append(row)
            b_vec.append(self.spoly.neg(equation.const))
        self._reduced_system = (tuple(tuple(row) for row in a_rows), tuple(b_vec), tuple(active_list), active_pos, substitutions)
        return self._reduced_system

    def reduceProbe(self, probe):
        _, _, active_list, active_pos, substitutions = self.prepareReducedSystem()
        out = probe.copy()
        for unknown_id, coeff in list(out.coeffs.items()):
            if unknown_id in active_pos:
                continue
            out.coeffs.pop(unknown_id, None)
            value = substitutions.get(unknown_id)
            if value is None:
                raise KernelSolveError("probe depends on an unknown removed from the reduced system without substitution")
            out = out.add(value.scale(coeff))
        for unknown_id in out.coeffs:
            if unknown_id not in active_pos:
                raise KernelSolveError("probe reduction left an inactive unknown")
        return out

    def probeNumerator(self, probe) -> SRat:
        a_rows, b_vec, active_list, active_pos, _ = self.prepareReducedSystem()
        reduced_probe = self.reduceProbe(probe)
        if not reduced_probe.coeffs and reduced_probe.const == self.spoly.zero:
            return self._zero_rat
        n = len(active_list)
        last_row = [self.spoly.zero for _ in range(n + 1)]
        for unknown_id, coeff in reduced_probe.coeffs.items():
            last_row[active_pos[unknown_id]] = self.spoly.neg(coeff)
        last_row[n] = reduced_probe.const
        matrix = [tuple(list(a_rows[row_id]) + [b_vec[row_id]]) for row_id in range(n)]
        matrix.append(tuple(last_row))
        return self.determinantRat(tuple(matrix))

    def determinantRat(self, matrix_key: tuple[tuple[int, ...], ...]) -> SRat:
        cached = self._det_cache.get(matrix_key)
        if cached is not None:
            return cached
        n = len(matrix_key)
        if n == 0:
            return self._one_rat
        if n <= 7:
            poly = self.determinantPolySparse(matrix_key)
            result = self.normalizeRat(SRat(poly, self.spoly.one))
            self._det_cache[matrix_key] = result
            return result
        matrix = [[self.ratFromPoly(value) for value in row] for row in matrix_key]
        sign = 1
        for k in range(n):
            pivot = self.choosePivot(matrix, k)
            if pivot is None:
                result = self._zero_rat
                self._det_cache[matrix_key] = result
                return result
            pivot_row, pivot_col = pivot
            if pivot_row != k:
                matrix[k], matrix[pivot_row] = matrix[pivot_row], matrix[k]
                sign *= -1
            if pivot_col != k:
                for row in matrix:
                    row[k], row[pivot_col] = row[pivot_col], row[k]
                sign *= -1
            pivot_value = matrix[k][k]
            for row_id in range(k + 1, n):
                entry = matrix[row_id][k]
                if entry.num == self.spoly.zero:
                    continue
                factor = self.ratDiv(entry, pivot_value)
                matrix[row_id][k] = self._zero_rat
                for col_id in range(k + 1, n):
                    matrix[row_id][col_id] = self.ratSub(matrix[row_id][col_id], self.ratMul(factor, matrix[k][col_id]))
        det = self._one_rat
        for idx in range(n):
            det = self.ratMul(det, matrix[idx][idx])
        if sign < 0:
            det = self.ratNeg(det)
        det = self.normalizeRat(det)
        self._det_cache[matrix_key] = det
        return det

    def determinantPolySparse(self, matrix_key: tuple[tuple[int, ...], ...]) -> int:
        n = len(matrix_key)
        nonzero_cols = {row_id: [col_id for col_id, value in enumerate(matrix_key[row_id]) if value != self.spoly.zero] for row_id in range(n)}
        assignment: dict[int, int] = {}

        def parity(cols: list[int]) -> int:
            value = 0
            for i in range(len(cols)):
                for j in range(i + 1, len(cols)):
                    if cols[i] > cols[j]:
                        value ^= 1
            return value

        def search(remaining_rows: set[int], used_cols: set[int], product: int) -> int:
            if not remaining_rows:
                cols = [assignment[row_id] for row_id in range(n)]
                return self.spoly.neg(product) if parity(cols) else product
            best_row = None
            best_cols: list[int] = []
            for row_id in remaining_rows:
                cols = [col_id for col_id in nonzero_cols[row_id] if col_id not in used_cols]
                if not cols:
                    return self.spoly.zero
                if best_row is None or len(cols) < len(best_cols):
                    best_row = row_id
                    best_cols = cols
            assert best_row is not None
            total = self.spoly.zero
            next_rows = set(remaining_rows)
            next_rows.remove(best_row)
            for col_id in best_cols:
                assignment[best_row] = col_id
                term = search(next_rows, used_cols | {col_id}, self.spoly.mul(product, matrix_key[best_row][col_id]))
                total = self.spoly.add(total, term)
                assignment.pop(best_row, None)
            return total

        return search(set(range(n)), set(), self.spoly.one)

    def evaluateRat(self, value: SRat, assignments: dict[str, Fraction], s_value: Fraction) -> Fraction:
        num = self.spoly.evaluate(value.num, assignments, s_value)
        den = self.spoly.evaluate(value.den, assignments, s_value)
        if den == 0:
            raise KernelError("numeric evaluation denominator is zero")
        return num / den

    def renderTex(self, response: Response, zero_pole: bool = False, fraction_symbol_threshold: int = 10) -> str:
        if response.value.den == self.spoly.zero:
            raise KernelZeroDenominatorError("cannot export response with zero denominator")
        if not isinstance(zero_pole, bool):
            raise ValueError("zero_pole must be True or False")
        if not isinstance(fraction_symbol_threshold, int) or fraction_symbol_threshold < 0:
            raise ValueError("fraction_symbol_threshold must be a nonnegative integer")
        value = self.normalizeForTex(self.applyPassiveValueAliases(response.value, response.passive_value_aliases))
        renderer = TexRenderer(self, fraction_symbol_threshold)
        rendered_value = renderer.preconditionRatForCoefficientRendering(value)
        renderer.prepareAggregateDefinitions(rendered_value)
        body, aggregate_definitions, coefficient_definitions = renderer.renderRatWithStableAggregates(rendered_value)
        lines = self.displayEquation("H(s)", body)
        if aggregate_definitions:
            lines.append("")
            for name, text in aggregate_definitions:
                lines.extend(self.displayEquation(name, text))
        if coefficient_definitions:
            lines.append("")
            for name, text in coefficient_definitions:
                lines.extend(self.displayEquation(name, text))
        if zero_pole:
            lines.append("")
            lines.extend(self.renderZeroPole(rendered_value))
        return "\n".join(lines) + "\n"

    def applyPassiveValueAliases(self, value: SRat, aliases: dict[str, str]) -> SRat:
        if not aliases or value.num == self.spoly.zero:
            return value
        cache: dict[int, int] = {}
        num = self.substitutePassiveValueAliasesInPoly(value.num, aliases, cache)
        den = self.substitutePassiveValueAliasesInPoly(value.den, aliases, cache)
        return self.normalizeRat(SRat(num, den))

    def substitutePassiveValueAliasesInPoly(self, poly_id: int, aliases: dict[str, str], cache: dict[int, int]) -> int:
        return self.spoly.poly(self.substitutePassiveValueAliasesInExpr(coeff, aliases, cache) for coeff in self.spoly.coeffs[poly_id])

    def substitutePassiveValueAliasesInExpr(self, expr_id: int, aliases: dict[str, str], cache: dict[int, int]) -> int:
        cached = cache.get(expr_id)
        if cached is not None:
            return cached

        kind = self.expr.kind[expr_id]
        if kind in {ZERO, ONE, RATIONAL}:
            result = expr_id
        elif kind == ATOM:
            atom_id = int(self.expr.value[expr_id])
            atom_kind = self.atoms.kind[atom_id]
            if atom_kind == "reciprocal_atom":
                base_id = self.atoms.base_param[atom_id]
                if base_id is None:
                    raise KernelError("reciprocal atom has no base parameter")
                alias = aliases.get(self.atoms.names[base_id])
                result = self.expr.atom(self.atoms.invParam(alias)) if alias else expr_id
            else:
                alias = aliases.get(self.atoms.names[atom_id])
                result = self.expr.atom(self.atoms.param(alias)) if alias else expr_id
        elif kind == ADD:
            result = self.expr.add(*(self.substitutePassiveValueAliasesInExpr(child, aliases, cache) for child in self.expr.nodeChildren(expr_id)))
        elif kind == MUL:
            result = self.expr.mul(*(self.substitutePassiveValueAliasesInExpr(child, aliases, cache) for child in self.expr.nodeChildren(expr_id)))
        elif kind == POW_INT:
            base, exp = self.expr.value[expr_id]
            result = self.expr.powInt(self.substitutePassiveValueAliasesInExpr(int(base), aliases, cache), int(exp))
        else:
            result = expr_id

        cache[expr_id] = result
        return result

    def displayEquation(self, lhs: str, rhs: str) -> list[str]:
        single_line = f"{lhs}={rhs}"
        if len(single_line) <= 5000:
            return [r"\[", single_line, r"\]"]
        rows = self.equationRows(lhs, rhs)
        return [r"\[", r"\begin{aligned}", *rows, r"\end{aligned}", r"\]"]

    def displayText(self, text: str) -> list[str]:
        return [r"\[", r"\text{" + text + r"}", r"\]"]

    def equationRows(self, lhs: str, rhs: str) -> list[str]:
        terms = self.splitTopLevelAddTerms(rhs)
        if not terms:
            return [f"{lhs}=0"]
        max_body_len = 5000
        groups: list[str] = []
        current = ""
        for term in terms:
            candidate = current + term if current else term
            if current and len(candidate) > max_body_len:
                groups.append(current)
                current = term
            else:
                current = candidate
        if current:
            groups.append(current)
        rows: list[str] = []
        for idx, group in enumerate(groups):
            if idx == 0:
                body = group[1:] if group.startswith("+") else group
                rows.append(f"{lhs}&={body}")
            else:
                rows.append(r"&{}" + group)
        return [row + (r"\\" if idx + 1 < len(rows) else "") for idx, row in enumerate(rows)]


    def splitTopLevelAddTerms(self, text: str) -> list[str]:
        terms: list[str] = []
        start = 0
        brace_depth = 0
        left_depth = 0
        idx = 0
        while idx < len(text):
            if text.startswith(r"\left", idx):
                left_depth += 1
                idx += 5
                continue
            if text.startswith(r"\right)", idx):
                left_depth = max(0, left_depth - 1)
                idx += 6
                continue
            char = text[idx]
            if char == "{":
                brace_depth += 1
            elif char == "}":
                brace_depth = max(0, brace_depth - 1)
            elif char in "+-" and idx > start and brace_depth == 0 and left_depth == 0:
                terms.append(text[start:idx])
                start = idx
            idx += 1
        terms.append(text[start:])
        return [term for term in terms if term]

    def normalizeForTex(self, value: SRat) -> SRat:
        if value.den == self.spoly.zero or value.num == self.spoly.zero:
            return value
        lead = self.spoly.leadingCoeff(value.den)
        if self.expr.isNegative(lead):
            return SRat(self.spoly.neg(value.num), self.spoly.neg(value.den))
        return value

    def renderZeroPole(self, value: SRat) -> list[str]:
        renderer = TexRenderer(self, 0)
        lines: list[str] = self.displayText("zero-pole analysis") + [""]
        zero_lines = self.rootDefinitions(value.num, "Z", renderer)
        pole_lines = self.rootDefinitions(value.den, "P", renderer)
        lines.extend(zero_lines or self.displayText("zeros: none"))
        lines.append("")
        lines.extend(pole_lines or self.displayText("poles: none"))
        return lines

    def rootDefinitions(self, poly_id: int, prefix: str, renderer: "TexRenderer") -> list[str]:
        extracted, residuals = self.extractPoleZeroFactors(poly_id)
        lines: list[str] = []
        index = 1
        for factor_id in extracted:
            degree = self.spoly.degree(factor_id)
            if degree == 1:
                root_text = self.linearRootText(factor_id, renderer)
                lines.extend(self.displayEquation(f"{prefix}_{{{index}}}", "s-" + root_text))
                index += 1
            elif degree == 2:
                first, second = self.quadraticRootTexts(factor_id, renderer)
                lines.extend(self.displayEquation(f"{prefix}_{{{index}}}", "s-" + first))
                index += 1
                lines.extend(self.displayEquation(f"{prefix}_{{{index}}}", "s-" + second))
                index += 1
            else:
                raise KernelTexError("internal pole-zero factor degree is not supported")
        if residuals:
            has_unresolved = any(self.spoly.degree(factor_id) > 0 for factor_id in residuals)
            if has_unresolved:
                label = "zeros" if prefix == "Z" else "poles"
                lines.extend(self.displayText(f"remaining {label}: analytic extraction failed"))
        return lines

    def extractPoleZeroFactors(self, poly_id: int) -> tuple[list[int], list[int]]:
        if self.spoly.degree(poly_id) <= 0:
            return [], []
        pending = [poly_id]
        extracted: list[int] = []
        residuals: list[int] = []
        s_factor = self.spoly.poly((self.expr.zero, self.expr.one))
        while pending:
            current = pending.pop()
            current = self.removeLeadingScalarOnly(current)
            if self.spoly.degree(current) <= 0:
                continue
            while self.spoly.degree(current) > 0 and self.spoly.coeffs[current][0] == self.expr.zero:
                extracted.append(s_factor)
                current = self.spoly.poly(self.spoly.coeffs[current][1:])
                current = self.removeLeadingScalarOnly(current)
            if self.spoly.degree(current) <= 0:
                continue
            known = self.spoly.known_factor_children.get(current)
            if known is not None and len(known) > 1:
                pending.extend(reversed(list(known)))
                continue
            degree = self.spoly.degree(current)
            if degree in (1, 2):
                extracted.append(current)
            else:
                residuals.append(current)
        return extracted, residuals

    def removeLeadingScalarOnly(self, poly_id: int) -> int:
        if poly_id == self.spoly.zero:
            return poly_id
        coeffs = self.spoly.coeffs[poly_id]
        if not coeffs:
            return poly_id
        lead = coeffs[-1]
        if lead == self.expr.one:
            return poly_id
        divided, ok = self.spoly.exactDiv(poly_id, self.spoly.poly((lead,)))
        return divided if ok else poly_id

    def linearRootText(self, factor_id: int, renderer: "TexRenderer") -> str:
        coeffs = self.spoly.coeffs[factor_id]
        a0 = coeffs[0] if len(coeffs) > 0 else self.expr.zero
        a1 = coeffs[1]
        if a0 == self.expr.zero:
            return "0"
        numerator = renderer.renderExprOptimized(self.expr.neg(a0))
        denominator = renderer.renderExprOptimized(a1)
        if denominator == "1":
            return self.parenthesizeRoot(numerator)
        return self.parenthesizeRoot(r"\frac{" + numerator + "}{" + denominator + "}")

    def quadraticRootTexts(self, factor_id: int, renderer: "TexRenderer") -> tuple[str, str]:
        coeffs = self.spoly.coeffs[factor_id]
        a0 = coeffs[0] if len(coeffs) > 0 else self.expr.zero
        a1 = coeffs[1] if len(coeffs) > 1 else self.expr.zero
        a2 = coeffs[2]
        a0_text = renderer.renderExprOptimized(a0)
        a1_text = renderer.renderExprOptimized(a1)
        a2_text = renderer.renderExprOptimized(a2)
        disc = a1_text + r"^{2}-4\left(" + a2_text + r"\right)\left(" + a0_text + r"\right)"
        root_plus = r"\frac{-\left(" + a1_text + r"\right)+\sqrt{" + disc + r"}}{2\left(" + a2_text + r"\right)}"
        root_minus = r"\frac{-\left(" + a1_text + r"\right)-\sqrt{" + disc + r"}}{2\left(" + a2_text + r"\right)}"
        return self.parenthesizeRoot(root_plus), self.parenthesizeRoot(root_minus)

    def parenthesizeRoot(self, text: str) -> str:
        if text == "0":
            return "0"
        return r"\left(" + text + r"\right)"


class TexRenderer:
    def __init__(self, kernel: SolverKernel, threshold: int):
        self.kernel = kernel
        self.expr = kernel.expr
        self.spoly = kernel.spoly
        self.threshold = threshold
        self.definitions: list[tuple[str, str]] = []
        self._definition_by_text: dict[str, str] = {}
        self._expr_cache: dict[int, str] = {}
        self._optimized_expr_cache: dict[int, str] = {}
        self._contains_rec_cache: dict[tuple[int, int], bool] = {}
        self._clear_rec_cache: dict[tuple[int, int], int] = {}
        self._expr_cost_cache: dict[int, int] = {}
        self._expr_symbol_cost_cache: dict[int, int] = {}
        self._coefficient_factor_cache: dict[int, int] = {}
        self._coefficient_normalize_cache: dict[int, int] = {}
        self._contains_expr_cache: dict[tuple[int, int], bool] = {}
        self._collect_linear_cache: dict[tuple[int, int], tuple[int, int]] = {}
        self._repeated_add_cache: dict[int, list[int]] = {}
        self._expand_cache: dict[tuple[int, int], int] = {}
        self._expr_metric_cache: dict[int, tuple[int, int, int]] = {}
        self._coefficient_expand_cache: dict[int, int] = {}
        self._coefficient_post_cache: dict[int, int] = {}
        self._aggregate_replace_cache: dict[int, int] = {}
        self._rational_absorb_cache: dict[int, int] = {}
        self._cancellable_reciprocal_cache: dict[int, int] = {}
        self._definition_scale_cache: dict[int, list[tuple[int, int]]] = {}
        self.aggregate_name_to_virtual_expr: dict[str, int] = {}
        self.aggregate_expr_to_name: dict[int, str] = {}
        self.aggregate_name_to_expr: dict[str, int] = {}
        self.aggregate_name_to_rhs: dict[str, str] = {}
        self.aggregate_name_to_atom_count: dict[str, int] = {}
        self.aggregate_name_to_atom_exprs: dict[str, tuple[int, ...]] = {}
        self.aggregate_name_to_kind: dict[str, str] = {}
        self.aggregate_disabled: set[str] = set()
        self._aggregate_use_counts: dict[str, int] = {}

    def preconditionRatForCoefficientRendering(self, value: SRat) -> SRat:
        """Choose a display-equivalent rational form before coefficient E_i creation.

        The solver keeps conductance-first atoms internally.  TeX rendering can
        still be shorter after multiplying the whole numerator and denominator
        by selected reciprocal bases and then dividing common display factors.
        This is renderer-local so solver normalization and intermediate MNA
        expressions are not enlarged.
        """
        if value.num == self.spoly.zero or value.den == self.spoly.zero:
            return value
        if not self.shouldPreconditionRat(value.num, value.den):
            return value

        best = value
        best_score = self.preconditionRatScore(best.num, best.den)

        candidates: list[SRat] = []
        for num_id, den_id in self.ratReciprocalClearCandidatePolys(value.num, value.den):
            candidates.append(SRat(num_id, den_id))

        for candidate in candidates[:16]:
            candidate = self.divideCommonDisplayFactors(candidate)
            score = self.preconditionRatScore(candidate.num, candidate.den)
            if score < best_score:
                best = candidate
                best_score = score

        self.resetRenderState()
        return best

    def shouldPreconditionRat(self, num_id: int, den_id: int) -> bool:
        """Bound display preconditioning to expressions where trial rendering is cheap.

        Large multi-transistor responses are handled by aggregate and E_i
        definitions.  Rendering them once without those definitions merely to
        score a reciprocal-cleared candidate can dominate export time, so the
        whole-rational preconditioner is restricted to small and medium
        polynomials.
        """
        coeffs = [
            coeff
            for poly_id in (num_id, den_id)
            for coeff in self.spoly.coeffs[poly_id]
            if coeff != self.expr.zero
        ]
        if len(coeffs) > 16:
            return False

        total_symbols = 0
        total_nodes = 0
        for coeff in coeffs:
            total_symbols += self.exprSymbolCost(coeff)
            if total_symbols > 320:
                return False
            total_nodes += self.exprTreeNodeCount(coeff, 1800)
            if total_nodes > 1800:
                return False
        return True

    def divideCommonDisplayFactors(self, value: SRat) -> SRat:
        num_id = value.num
        den_id = value.den
        for factor in self.commonPolyFactorOrder(num_id, den_id):
            trial_num = self.divPolyByExpr(num_id, factor)
            trial_den = self.divPolyByExpr(den_id, factor)
            if trial_num is None or trial_den is None:
                continue
            num_id, den_id = trial_num, trial_den
        return SRat(num_id, den_id)

    def preconditionRatScore(self, num_id: int, den_id: int) -> tuple[int, int, int, int]:
        saved_threshold = self.threshold
        saved_definitions = list(self.definitions)
        saved_definition_by_text = dict(self._definition_by_text)
        self.threshold = 0
        try:
            text = self.renderRatFromPolys(num_id, den_id)
            return (
                self.fractionNestingPenalty(text) + self.longSimpleDenominatorPenalty(text),
                self.renderedSymbolCount(text),
                len(text),
                self.ratReciprocalOccurrence(num_id, den_id),
            )
        finally:
            self.threshold = saved_threshold
            self.definitions = saved_definitions
            self._definition_by_text = saved_definition_by_text

    def ratReciprocalClearCandidatePolys(self, num_id: int, den_id: int) -> list[tuple[int, int]]:
        atoms = self.polyReciprocalAtomExprs(num_id, den_id)
        if not atoms:
            return []

        limited = atoms[:8]
        specs: list[tuple[tuple[int, int], ...]] = []
        seen: set[tuple[tuple[int, int], ...]] = set()

        def addSpec(items: list[tuple[int, int]]) -> None:
            spec = tuple((atom, count) for atom, count in items if count > 0)
            if spec and spec not in seen:
                seen.add(spec)
                specs.append(spec)

        max_power_items = [
            (atom, self.polyReciprocalMaxPower(num_id, den_id, atom))
            for atom in limited
        ]
        addSpec(max_power_items)
        for atom, max_power in max_power_items:
            addSpec([(atom, max_power)])
        addSpec(max_power_items[:4])

        out: list[tuple[int, int]] = []
        for spec in specs[:12]:
            trial_num = num_id
            trial_den = den_id
            for atom_expr, count in spec:
                for _ in range(count):
                    trial_num = self.clearOneReciprocalInPoly(trial_num, atom_expr)
                    trial_den = self.clearOneReciprocalInPoly(trial_den, atom_expr)
            if trial_num != num_id or trial_den != den_id:
                out.append((trial_num, trial_den))
        return out

    def clearOneReciprocalInPoly(self, poly_id: int, atom_expr: int) -> int:
        return self.spoly.poly(
            self.clearOneReciprocal(coeff, atom_expr)
            for coeff in self.spoly.coeffs[poly_id]
        )

    def polyReciprocalAtomExprs(self, *poly_ids: int) -> list[int]:
        out: set[int] = set()
        for poly_id in poly_ids:
            for coeff in self.spoly.coeffs[poly_id]:
                self.collectReciprocalAtomExprs(coeff, out)
        return sorted(
            out,
            key=lambda atom_expr: -sum(
                self.reciprocalOccurrence(coeff, atom_expr)
                for poly_id in poly_ids
                for coeff in self.spoly.coeffs[poly_id]
            ),
        )

    def polyReciprocalMaxPower(self, num_id: int, den_id: int, atom_expr: int) -> int:
        return max(
            (
                self.reciprocalMaxPower(coeff, atom_expr)
                for poly_id in (num_id, den_id)
                for coeff in self.spoly.coeffs[poly_id]
            ),
            default=0,
        )

    def ratReciprocalOccurrence(self, num_id: int, den_id: int) -> int:
        atoms = self.polyReciprocalAtomExprs(num_id, den_id)
        return sum(
            self.reciprocalOccurrence(coeff, atom)
            for atom in atoms
            for poly_id in (num_id, den_id)
            for coeff in self.spoly.coeffs[poly_id]
        )

    def prepareAggregateDefinitions(self, value: SRat) -> None:
        self.aggregate_expr_to_name.clear()
        self.aggregate_name_to_expr.clear()
        self.aggregate_name_to_rhs.clear()
        self.aggregate_name_to_atom_count.clear()
        self.aggregate_name_to_atom_exprs.clear()
        self.aggregate_name_to_kind.clear()
        self.aggregate_name_to_virtual_expr.clear()
        self.aggregate_disabled.clear()

        all_cap_groups: dict[str, set[int]] = {}
        all_conductance_groups: dict[str, set[int]] = {}
        device_subset_counts: dict[tuple[str, str, tuple[int, ...]], int] = {}
        eff_subset_counts: dict[tuple[str, tuple[int, ...]], int] = {}
        for atom_expr in self.atomExprsInRat(value):
            kind_device = self.aggregateKindAndDevice(atom_expr)
            if kind_device is None:
                continue
            kind, device = kind_device
            target = all_cap_groups if kind == "C" else all_conductance_groups
            target.setdefault(device, set()).add(atom_expr)

        for poly_id in (value.num, value.den):
            for coeff in self.spoly.coeffs[poly_id]:
                if coeff != self.expr.zero:
                    self.collectAggregateSubsetCandidates(coeff, device_subset_counts, eff_subset_counts)

        def chooseDeviceGroups(prefix: str, all_groups: dict[str, set[int]]) -> dict[str, tuple[int, ...]]:
            chosen: dict[str, tuple[int, ...]] = {}
            for device, atom_set in all_groups.items():
                candidates: dict[tuple[int, ...], int] = {}
                full = tuple(sorted(atom_set, key=self.expr.sortKey))
                if len(full) >= 2:
                    candidates[full] = candidates.get(full, 0) + 1
                for (kind, cand_device, subset), count in device_subset_counts.items():
                    if kind == prefix and cand_device == device and len(subset) >= 2:
                        if set(subset).issubset(atom_set):
                            candidates[subset] = candidates.get(subset, 0) + count
                if not candidates:
                    continue
                def score(item: tuple[tuple[int, ...], int]) -> tuple[int, int, int]:
                    subset, count = item
                    k = len(subset)
                    return (count * (k - 1) - k, count * (k - 1), k)
                best_subset, best_count = max(candidates.items(), key=lambda item: (score(item), tuple(self.expr.sortKey(x) for x in item[0])))
                if len(best_subset) >= 2 and best_count > 0:
                    chosen[device] = best_subset
            return chosen

        cap_groups = chooseDeviceGroups("C", all_cap_groups)
        conductance_groups = chooseDeviceGroups("G", all_conductance_groups)
        for prefix, groups in (("C", cap_groups), ("G", conductance_groups)):
            for device in sorted(groups):
                atom_exprs = tuple(sorted(groups[device], key=self.expr.sortKey))
                if len(atom_exprs) < 2:
                    continue
                self.addAggregateDefinition(prefix, f"{prefix}_{{{device}}}", atom_exprs)

        used_eff_subsets = self.chooseEffectiveAggregateSubsets(eff_subset_counts)
        counters = {"C": 0, "G": 0, "R": 0, "L": 0}
        for kind, atom_exprs in used_eff_subsets:
            counters[kind] += 1
            name = f"{kind}_{{eff{counters[kind]}}}"
            self.addAggregateDefinition(kind + "eff", name, atom_exprs)

    def addAggregateDefinition(self, kind: str, name: str, atom_exprs: tuple[int, ...]) -> None:
        if len(atom_exprs) < 2:
            return
        atom_exprs = tuple(sorted(atom_exprs, key=self.expr.sortKey))
        expr_id = self.expr.add(*atom_exprs)
        if expr_id in self.aggregate_expr_to_name:
            return
        rhs = self.renderAggregateDefinitionRhs(kind, atom_exprs)
        self.aggregate_expr_to_name[expr_id] = name
        self.aggregate_name_to_expr[name] = expr_id
        self.aggregate_name_to_rhs[name] = rhs
        self.aggregate_name_to_atom_count[name] = len(atom_exprs)
        self.aggregate_name_to_atom_exprs[name] = atom_exprs
        self.aggregate_name_to_kind[name] = kind

    def chooseEffectiveAggregateSubsets(
        self,
        subset_counts: dict[tuple[str, tuple[int, ...]], int],
    ) -> list[tuple[str, tuple[int, ...]]]:
        candidates: list[tuple[int, int, str, tuple[int, ...]]] = []
        existing_exprs = set(self.aggregate_expr_to_name)
        for (kind, subset), count in subset_counts.items():
            if len(subset) < 2:
                continue
            expr_id = self.expr.add(*subset)
            if expr_id in existing_exprs:
                continue
            k = len(subset)
            score = count * (k - 1) - (k + 1)
            if score <= 0:
                continue
            candidates.append((score, count, kind, subset))
        candidates.sort(key=lambda item: (-item[0], -item[1], item[2], tuple(self.expr.sortKey(x) for x in item[3])))
        selected: list[tuple[str, tuple[int, ...]]] = []
        selected_exprs: set[int] = set(existing_exprs)
        per_kind_count = {"C": 0, "G": 0, "R": 0, "L": 0}
        for score, count, kind, subset in candidates:
            if per_kind_count[kind] >= 24:
                continue
            expr_id = self.expr.add(*subset)
            if expr_id in selected_exprs:
                continue
            selected.append((kind, subset))
            selected_exprs.add(expr_id)
            per_kind_count[kind] += 1
            if len(selected) >= 64:
                break
        return selected

    def collectAggregateSubsetCandidates(
        self,
        expr_id: int,
        device_subset_counts: dict[tuple[str, str, tuple[int, ...]], int],
        eff_subset_counts: dict[tuple[str, tuple[int, ...]], int],
        visited: set[int] | None = None,
    ) -> None:
        if visited is None:
            visited = set()
        if expr_id in visited:
            return
        visited.add(expr_id)
        kind = self.expr.kind[expr_id]
        if kind == ADD:
            by_kind_device_coeff: dict[tuple[str, str, Fraction], set[int]] = {}
            by_eff_kind_coeff: dict[tuple[str, Fraction], set[int]] = {}
            for child in self.expr.nodeChildren(expr_id):
                coeff, base = self.expr.splitCoeff(child)
                kind_device = self.aggregateKindAndDevice(base)
                if kind_device is not None:
                    cap_or_g, device = kind_device
                    by_kind_device_coeff.setdefault((cap_or_g, device, coeff), set()).add(base)
                eff_kind = self.effectiveAggregateKind(base)
                if eff_kind is not None:
                    by_eff_kind_coeff.setdefault((eff_kind, coeff), set()).add(base)
            for (cap_or_g, device, _coeff), atom_exprs in by_kind_device_coeff.items():
                if len(atom_exprs) >= 2:
                    subset = tuple(sorted(atom_exprs, key=self.expr.sortKey))
                    key = (cap_or_g, device, subset)
                    device_subset_counts[key] = device_subset_counts.get(key, 0) + 1
            for (eff_kind, _coeff), atom_exprs in by_eff_kind_coeff.items():
                if len(atom_exprs) >= 2:
                    subset = tuple(sorted(atom_exprs, key=self.expr.sortKey))
                    key = (eff_kind, subset)
                    eff_subset_counts[key] = eff_subset_counts.get(key, 0) + 1
        for child in self.expr.nodeChildren(expr_id):
            self.collectAggregateSubsetCandidates(child, device_subset_counts, eff_subset_counts, visited)

    def atomExprsInRat(self, value: SRat) -> set[int]:
        out: set[int] = set()
        for poly_id in (value.num, value.den):
            for coeff in self.spoly.coeffs[poly_id]:
                self.collectAtomExprs(coeff, out)
        return out

    def collectAtomExprs(self, expr_id: int, out: set[int]) -> None:
        kind = self.expr.kind[expr_id]
        if kind == ATOM:
            out.add(expr_id)
            return
        if kind in {ZERO, ONE, RATIONAL}:
            return
        for child in self.expr.nodeChildren(expr_id):
            self.collectAtomExprs(child, out)

    def aggregateKindAndDevice(self, atom_expr: int) -> tuple[str, str] | None:
        if self.expr.kind[atom_expr] != ATOM:
            return None
        atom_id = int(self.expr.value[atom_expr])
        kind = self.expr.atoms.kind[atom_id]
        name = self.expr.atoms.names[atom_id]
        if kind == "reciprocal_atom":
            base_id = self.expr.atoms.base_param[atom_id]
            if base_id is None:
                return None
            base_name = self.expr.atoms.names[base_id]
            if not base_name.startswith("r_"):
                return None
            device = self.deviceNameFromParam(base_name)
            return ("G", device) if self.isActiveDeviceName(device) else None
        if name.startswith("C_"):
            device = self.deviceNameFromParam(name)
            return ("C", device) if self.isActiveDeviceName(device) else None
        if name.startswith("g_"):
            device = self.deviceNameFromParam(name)
            return ("G", device) if self.isActiveDeviceName(device) else None
        return None

    def effectiveAggregateKind(self, atom_expr: int) -> str | None:
        if self.expr.kind[atom_expr] != ATOM:
            return None
        atom_id = int(self.expr.value[atom_expr])
        kind = self.expr.atoms.kind[atom_id]
        name = self.expr.atoms.names[atom_id]
        if kind == "reciprocal_atom":
            base_id = self.expr.atoms.base_param[atom_id]
            if base_id is None:
                return None
            base_name = self.expr.atoms.names[base_id]
            if base_name.startswith(("r_", "R_")):
                return "G"
            if base_name.startswith("g_"):
                return "R"
            return None
        if name.startswith("C_"):
            return "C"
        if name.startswith("g_"):
            return "G"
        if name.startswith(("r_", "R_")):
            return "R"
        if name.startswith("L_"):
            return "L"
        return None

    def deviceNameFromParam(self, name: str) -> str:
        rest = name.split("_", 1)[1] if "_" in name else name
        return rest.split("_", 1)[0]

    def isActiveDeviceName(self, device: str) -> bool:
        low = device.lower()
        return low.startswith(("mn", "mp", "m", "qn", "qp", "q"))

    def renderAggregateDefinitionRhs(self, kind: str, atom_exprs: tuple[int, ...]) -> str:
        raw = self.joinAggregateDefinitionTerms(atom_exprs)
        if kind not in {"Ceff", "Geff"}:
            return raw
        base_kind = "C" if kind == "Ceff" else "G"
        remaining = set(atom_exprs)
        chosen: list[str] = []
        candidates: list[tuple[int, str, tuple[int, ...]]] = []
        for name, exprs in self.aggregate_name_to_atom_exprs.items():
            if self.aggregate_name_to_kind.get(name) != base_kind:
                continue
            expr_set = set(exprs)
            if len(exprs) >= 2 and expr_set.issubset(remaining):
                candidates.append((len(exprs) - 1, name, exprs))
        candidates.sort(key=lambda item: (-item[0], item[1]))
        for _, name, exprs in candidates:
            expr_set = set(exprs)
            if expr_set.issubset(remaining):
                chosen.append(name)
                remaining.difference_update(expr_set)
        if not chosen:
            return raw
        sorted_remaining = sorted(remaining, key=self.expr.sortKey)
        grouped_texts = chosen + [self.renderExprNoAggregates(atom_expr) for atom_expr in sorted_remaining]
        grouped = self.joinSignedTerms(grouped_texts)
        return grouped if self.renderedSymbolCount(grouped) < self.renderedSymbolCount(raw) else raw

    def joinAggregateDefinitionTerms(self, atom_exprs: tuple[int, ...]) -> str:
        parts: list[str] = []
        for atom_expr in atom_exprs:
            text = self.renderExprNoAggregates(atom_expr)
            if text.startswith("-"):
                parts.append(text)
            elif parts:
                parts.append("+" + text)
            else:
                parts.append(text)
        return "".join(parts)

    def renderExprNoAggregates(self, expr_id: int) -> str:
        disabled = self.aggregate_disabled
        saved = set(disabled)
        self.aggregate_disabled = set(self.aggregate_name_to_expr)
        try:
            self._expr_cache.clear()
            self._optimized_expr_cache.clear()
            return self.renderExpr(expr_id)
        finally:
            self.aggregate_disabled = saved
            self._expr_cache.clear()
            self._optimized_expr_cache.clear()

    def renderRatWithStableAggregates(self, value: SRat) -> tuple[str, list[tuple[str, str]], list[tuple[str, str]]]:
        previous_disabled: set[str] | None = None
        body = ""
        coefficient_definitions: list[tuple[str, str]] = []
        used_aggregates: set[str] = set()
        for _ in range(8):
            self.resetRenderState()
            body = self.renderRat(value)
            coefficient_definitions = list(self.definitions)
            all_text = body + "\n" + "\n".join(rhs for _, rhs in coefficient_definitions)
            directly_used = self.aggregateNamesUsedInText(all_text)
            used_aggregates = self.aggregateDependencyClosure(directly_used)
            disabled = set(self.aggregate_disabled)
            for name in sorted(directly_used):
                if not self.aggregateUseShortensTotal(name, all_text):
                    disabled.add(name)
            if disabled == self.aggregate_disabled or disabled == previous_disabled:
                break
            previous_disabled = set(self.aggregate_disabled)
            self.aggregate_disabled = disabled
        aggregate_definitions = [(name, self.aggregate_name_to_rhs[name]) for name in sorted(used_aggregates, key=self.aggregateSortKey) if name not in self.aggregate_disabled]
        return body, aggregate_definitions, coefficient_definitions

    def resetRenderState(self) -> None:
        self.definitions.clear()
        self._definition_by_text.clear()
        self._expr_cache.clear()
        self._optimized_expr_cache.clear()
        self._contains_rec_cache.clear()
        self._clear_rec_cache.clear()
        self._expr_cost_cache.clear()
        self._expr_symbol_cost_cache.clear()
        self._coefficient_factor_cache.clear()
        self._coefficient_normalize_cache.clear()
        self._contains_expr_cache.clear()
        self._collect_linear_cache.clear()
        self._repeated_add_cache.clear()
        self._coefficient_expand_cache.clear()
        self._coefficient_post_cache.clear()
        self._aggregate_replace_cache.clear()
        self._aggregate_use_counts.clear()
        self._definition_scale_cache.clear()

    def aggregateNamesUsedInText(self, text: str) -> set[str]:
        out: set[str] = set()
        for name in self.aggregate_name_to_expr:
            if name not in self.aggregate_disabled and name in text:
                out.add(name)
        return out

    def aggregateDependencyClosure(self, names: set[str]) -> set[str]:
        out = set(names)
        changed = True
        while changed:
            changed = False
            for name in list(out):
                rhs = self.aggregate_name_to_rhs.get(name, "")
                for dep in self.aggregateNamesUsedInText(rhs):
                    if dep not in out:
                        out.add(dep)
                        changed = True
        return out

    def aggregateUseShortensTotal(self, name: str, text: str) -> bool:
        uses = text.count(name)
        if uses <= 0:
            return False
        atom_count = self.aggregate_name_to_atom_count[name]
        # Each circuit symbol has unit length. A used aggregate occurrence costs one symbol;
        # the definition costs one LHS symbol plus the RHS atom symbols.
        return uses * atom_count > uses + atom_count + 1

    def aggregateSortKey(self, name: str) -> tuple[int, int, str]:
        kind = self.aggregate_name_to_kind.get(name, "")
        order = {"C": 0, "G": 1, "Ceff": 2, "Geff": 3, "Reff": 4, "Leff": 5}.get(kind, 9)
        match = re.search(r"eff(\d+)", name)
        suffix = int(match.group(1)) if match else 0
        return (order, suffix, name)

    def aggregateLinearMatch(self, expr_id: int) -> tuple[str, Fraction] | None:
        exact_name = self.aggregate_expr_to_name.get(expr_id)
        if exact_name is not None and exact_name not in self.aggregate_disabled:
            return exact_name, Fraction(1)
        coeff, base = self.expr.splitCoeff(expr_id)
        if base != expr_id:
            match = self.aggregateLinearMatch(base)
            if match is not None:
                name, base_coeff = match
                return name, coeff * base_coeff
        if self.expr.kind[expr_id] != ADD:
            return None
        by_base: dict[int, Fraction] = {}
        for child in self.expr.nodeChildren(expr_id):
            child_coeff, child_base = self.expr.splitCoeff(child)
            if child_base == self.expr.one:
                return None
            by_base[child_base] = by_base.get(child_base, Fraction(0)) + child_coeff
        if not by_base:
            return None
        coeffs = {value for value in by_base.values() if value != 0}
        if len(coeffs) != 1:
            return None
        common_coeff = next(iter(coeffs))
        base_set = {base for base, value in by_base.items() if value != 0}
        for name, atom_exprs in self.aggregate_name_to_atom_exprs.items():
            if name in self.aggregate_disabled:
                continue
            if base_set == set(atom_exprs):
                return name, common_coeff
        return None

    def isAggregateLinearExpr(self, expr_id: int) -> bool:
        return self.aggregateLinearMatch(expr_id) is not None

    def aggregateVirtualExpr(self, name: str) -> int:
        existing = self.aggregate_name_to_virtual_expr.get(name)
        if existing is not None:
            return existing
        atom_id = self.expr.atoms.param(name)
        expr_id = self.expr.atom(atom_id)
        self.aggregate_name_to_virtual_expr[name] = expr_id
        return expr_id

    def replaceAggregateLinearExprs(self, expr_id: int) -> int:
        cached = self._aggregate_replace_cache.get(expr_id)
        if cached is not None:
            return cached
        match = self.aggregateLinearMatch(expr_id)
        if match is not None:
            name, coeff = match
            result = self.expr.mul(self.expr.rational(coeff), self.aggregateVirtualExpr(name))
            self._aggregate_replace_cache[expr_id] = result
            return result
        kind = self.expr.kind[expr_id]
        if kind in {ZERO, ONE, RATIONAL, ATOM}:
            result = expr_id
        elif kind == ADD:
            result = self.expr.add(*(self.replaceAggregateLinearExprs(child) for child in self.expr.nodeChildren(expr_id)))
        elif kind == MUL:
            result = self.expr.mul(*(self.replaceAggregateLinearExprs(child) for child in self.expr.nodeChildren(expr_id)))
        elif kind == POW_INT:
            base, exp = self.expr.value[expr_id]
            result = self.expr.powInt(self.replaceAggregateLinearExprs(int(base)), int(exp))
        else:
            result = expr_id
        self._aggregate_replace_cache[expr_id] = result
        return result

    def renderRat(self, value: SRat) -> str:
        best = self.renderRatFromPolys(value.num, value.den)
        num_id, den_id = value.num, value.den
        for factor in self.commonPolyFactorOrder(num_id, den_id):
            trial_num = self.divPolyByExpr(num_id, factor)
            trial_den = self.divPolyByExpr(den_id, factor)
            if trial_num is None or trial_den is None:
                continue
            candidate = self.renderRatFromPolys(trial_num, trial_den)
            if self.betterText(candidate, best):
                best = candidate
                num_id, den_id = trial_num, trial_den
        return best

    def renderRatFromPolys(self, num_id: int, den_id: int) -> str:
        num = self.renderPoly(num_id)
        den = self.renderPoly(den_id)
        if den_id == self.spoly.one:
            return num
        return self.fractionText(num, den)

    def renderPoly(self, poly_id: int) -> str:
        raw = self.renderPolyExpanded(poly_id)
        best = raw
        for factor in self.polyCommonFactorOrder(poly_id):
            divided = self.divPolyByExpr(poly_id, factor)
            if divided is None or divided == poly_id:
                continue
            factor_text = self.renderExprOptimized(factor)
            divided_text = self.renderPolyExpanded(divided)
            candidate = self.productText(factor_text, divided_text)
            if self.betterText(candidate, best):
                best = candidate
        return best

    def renderPolyExpanded(self, poly_id: int) -> str:
        if poly_id == self.spoly.zero:
            return "0"
        coeffs = self.spoly.coeffs[poly_id]
        terms: list[str] = []
        for power in range(len(coeffs) - 1, -1, -1):
            coeff = coeffs[power]
            if coeff == self.expr.zero:
                continue
            coeff_text = self.renderCoeff(coeff)
            negative = coeff_text.startswith("-")
            coeff_body = coeff_text[1:] if negative else coeff_text
            term = self.termText(coeff_body, power)
            if not terms:
                terms.append(("-" if negative else "") + term)
            else:
                terms.append(("-" if negative else "+") + term)
        return "".join(terms) if terms else "0"

    def divPolyByExpr(self, poly_id: int, divisor: int) -> int | None:
        out: list[int] = []
        for coeff in self.spoly.coeffs[poly_id]:
            quotient = self.expr.divByExpr(coeff, divisor)
            if quotient is None:
                return None
            out.append(quotient)
        return self.spoly.poly(out)

    def polyCommonFactorMap(self, poly_id: int) -> dict[int, int]:
        common: dict[int, int] | None = None
        for coeff in self.spoly.coeffs[poly_id]:
            if coeff == self.expr.zero:
                continue
            current_items = self.expr.commonTermFactorPowerItems(coeff)
            if common is None:
                common = dict(current_items)
                continue
            current = dict(current_items)
            for factor in list(common):
                if factor not in current:
                    del common[factor]
                else:
                    common[factor] = min(common[factor], current[factor])
            if not common:
                break
        return common if common is not None else {}

    def polyCommonFactorOrder(self, poly_id: int) -> list[int]:
        factors = self.polyCommonFactorMap(poly_id)
        candidates: list[int] = []
        for factor, power in factors.items():
            if factor in (self.expr.zero, self.expr.one) or power <= 0:
                continue
            candidates.append(self.expr.powInt(factor, power) if power > 1 else factor)
        candidates.sort(key=lambda item: (-len(self.renderExpr(item)), item))
        if len(candidates) > 1:
            product = self.expr.mul(*candidates)
            candidates.insert(0, product)
        return candidates[:12]

    def commonPolyFactorOrder(self, num_id: int, den_id: int) -> list[int]:
        left = self.polyCommonFactorMap(num_id)
        right = self.polyCommonFactorMap(den_id)
        common: dict[int, int] = {}
        for factor, power in left.items():
            if factor in right:
                common[factor] = min(power, right[factor])
        candidates: list[int] = []
        for factor, power in common.items():
            if factor in (self.expr.zero, self.expr.one) or power <= 0:
                continue
            for _ in range(power):
                candidates.append(factor)
        candidates.sort(key=lambda item: (-len(self.renderExpr(item)), item))
        return candidates[:32]

    def productText(self, left: str, right: str) -> str:
        if left == "1":
            return right
        if right == "1":
            return left
        if self.needsParens(right):
            right = r"\left(" + right + r"\right)"
        if self.needsParens(left):
            left = r"\left(" + left + r"\right)"
        return left + " " + right

    def termText(self, coeff: str, power: int) -> str:
        if power == 0:
            return coeff
        s_part = "s" if power == 1 else f"s^{{{power}}}"
        if coeff == "1":
            return s_part
        if self.needsParens(coeff):
            return s_part + r"\left(" + coeff + r"\right)"
        return s_part + " " + coeff

    def renderCoeff(self, expr_id: int) -> str:
        coeff_expr = self.normalizeCoefficientExpr(expr_id)
        coeff_expr = self.replaceAggregateLinearExprs(coeff_expr)
        coeff_expr = self.postOptimizeCoefficientExpr(coeff_expr)
        text = self.renderExprOptimized(coeff_expr)
        if self.threshold > 0 and self.renderedSymbolCount(text) > self.threshold:
            existing = self._definition_by_text.get(text)
            if existing is not None:
                return existing
            name = f"E_{{{len(self.definitions) + 1}}}"
            self._definition_by_text[text] = name
            lhs, rhs = self.definitionEquation(name, coeff_expr, text)
            self.definitions.append((lhs, rhs))
            return name
        return text

    def definitionEquation(self, name: str, expr_id: int, default_text: str) -> tuple[str, str]:
        if self.exprTreeNodeCount(expr_id, 2200) > 2200 or self.exprSymbolCost(expr_id) > 520:
            return name, default_text
        best_scale = self.expr.one
        best_rhs_expr = expr_id
        best_rhs_text = default_text
        best_score = self.definitionEquationScore(name, default_text)
        for scale_expr, rhs_expr in self.definitionScaleCandidates(expr_id):
            if scale_expr == self.expr.one:
                continue
            rhs_text = self.renderExprOptimized(rhs_expr)
            score = self.definitionEquationScore(self.definitionLhsText(name, scale_expr), rhs_text)
            if score < best_score:
                best_scale = scale_expr
                best_rhs_expr = rhs_expr
                best_rhs_text = rhs_text
                best_score = score
        if best_scale == self.expr.one:
            return name, best_rhs_text
        return name, self.renderScaledDefinitionRhs(best_scale, best_rhs_expr, best_rhs_text)

    def definitionEquationScore(self, lhs: str, rhs: str) -> tuple[int, int, int]:
        text = lhs + "=" + rhs
        return (
            self.renderedSymbolCount(text),
            self.fractionNestingPenalty(text) + self.longSimpleDenominatorPenalty(text),
            len(text),
        )

    def definitionLhsText(self, name: str, scale_expr: int) -> str:
        if scale_expr == self.expr.one:
            return name
        if self.expr.kind[scale_expr] == RATIONAL:
            p, q = self.expr.value[scale_expr]
            if int(q) == 1:
                value = int(p)
                if value == 1:
                    return name
                if value == -1:
                    return "-" + name
                return f"{value}{name}"
        scale_text = self.renderExprOptimized(scale_expr)
        if scale_text == "1":
            return name
        if scale_text == "-1":
            return "-" + name
        return self.productText(scale_text, name)

    def renderScaledDefinitionRhs(self, scale_expr: int, rhs_expr: int, rhs_text: str | None = None) -> str:
        rhs = rhs_text if rhs_text is not None else self.renderExprOptimized(rhs_expr)
        scale = self.renderExprOptimized(scale_expr)
        if scale == "1":
            return rhs
        if scale == "-1":
            return "-" + (rhs if not self.needsParens(rhs) else r"\left(" + rhs + r"\right)")
        body = rhs if not self.needsParens(rhs) else r"\left(" + rhs + r"\right)"
        if self.expr.kind[scale_expr] == RATIONAL:
            p, q = self.expr.value[scale_expr]
            if int(q) == 1 and int(p) > 0:
                return self.productText(self.fractionText("1", str(int(p))), body)
        return self.productText(self.fractionText("1", scale), body)

    def definitionScaleCandidates(self, expr_id: int) -> list[tuple[int, int]]:
        cached = self._definition_scale_cache.get(expr_id)
        if cached is not None:
            return cached
        out: list[tuple[int, int]] = []
        seen: set[tuple[int, int]] = set()

        def addCandidate(divisor: int, scale: int, numerator_factor: int) -> None:
            if divisor == self.expr.one or scale == self.expr.one:
                return
            quotient = self.expr.divByExpr(expr_id, divisor)
            if quotient is None or quotient == expr_id:
                return
            rhs = self.expr.mul(numerator_factor, quotient)
            key = (scale, rhs)
            if key not in seen:
                seen.add(key)
                out.append(key)

        coeff, base = self.expr.splitCoeff(expr_id)
        if coeff.denominator != 1:
            divisor = self.expr.rational(Fraction(1, coeff.denominator))
            scale = self.expr.rational(coeff.denominator)
            numerator_factor = self.expr.rational(coeff.numerator)
            addCandidate(divisor, scale, numerator_factor)
        if coeff.numerator != 0 and abs(coeff.numerator) != 1 and coeff.denominator != 1:
            divisor = self.expr.rational(coeff)
            scale = self.expr.rational(coeff.denominator)
            numerator_factor = self.expr.rational(coeff.numerator)
            addCandidate(divisor, scale, numerator_factor)

        rational_divisor = self.commonRationalDivisor(expr_id)
        if rational_divisor is not None:
            p, q = self.expr.value[rational_divisor]
            addCandidate(rational_divisor, self.expr.rational(int(q)), self.expr.rational(int(p)))

        rec_factors: list[int] = []
        common_factors = self.expr.commonTermFactorPowerItems(expr_id)
        for factor, power in common_factors:
            if power <= 0 or not self.isReciprocalAtomExpr(factor):
                continue
            rec_factors.extend([factor] * min(power, 3))
        if rec_factors:
            rec_factors = rec_factors[:4]
            rec_divisor = self.expr.mul(*rec_factors)
            rec_scale = self.expr.mul(*(self.baseExprForReciprocalAtom(factor) for factor in rec_factors))
            addCandidate(rec_divisor, rec_scale, self.expr.one)
            if rational_divisor is not None:
                p, q = self.expr.value[rational_divisor]
                combo_divisor = self.expr.mul(rational_divisor, rec_divisor)
                combo_scale = self.expr.mul(self.expr.rational(int(q)), rec_scale)
                addCandidate(combo_divisor, combo_scale, self.expr.rational(int(p)))

        out.sort(key=lambda item: (self.exprSymbolCost(item[1]), self.exprTreeNodeCount(item[1], 100000), self.expr.sortKey(item[0])))
        self._definition_scale_cache[expr_id] = out[:8]
        return out

    def longSimpleDenominatorPenalty(self, text: str) -> int:
        penalty = 0
        start = 0
        needle = r"\frac{"
        while True:
            pos = text.find(needle, start)
            if pos < 0:
                break
            numerator_start = pos + len(needle)
            numerator_end = self.matchingBraceEnd(text, numerator_start - 1)
            if numerator_end < 0 or numerator_end + 1 >= len(text) or text[numerator_end + 1:numerator_end + 3] != "{":
                start = numerator_start
                continue
            denominator_start = numerator_end + 2
            denominator_end = self.matchingBraceEnd(text, denominator_start - 1)
            if denominator_end < 0:
                start = numerator_start
                continue
            numerator = text[numerator_start:numerator_end]
            denominator = text[denominator_start:denominator_end].strip()
            if denominator.isdigit() or re.fullmatch(r"[A-Za-z]+(?:_\{[^}]+\})?", denominator):
                symbol_count = self.renderedSymbolCount(numerator)
                if symbol_count >= 8:
                    penalty += symbol_count
            start = denominator_end + 1
        return penalty

    def matchingBraceEnd(self, text: str, open_index: int) -> int:
        if open_index < 0 or open_index >= len(text) or text[open_index] != "{":
            return -1
        depth = 0
        for index in range(open_index, len(text)):
            char = text[index]
            if char == "{":
                depth += 1
            elif char == "}":
                depth -= 1
                if depth == 0:
                    return index
        return -1

    def postOptimizeCoefficientExpr(self, expr_id: int) -> int:
        cached = self._coefficient_post_cache.get(expr_id)
        if cached is not None:
            return cached
        best_expr = expr_id
        best_text = self.renderExprOptimized(best_expr)
        candidates = [expr_id]

        # These passes are coefficient-local: after E_i has been chosen as an
        # s-polynomial coefficient, they only rewrite that coefficient.  They
        # first remove algebraically cancellable rational/reciprocal clutter,
        # then search common factors on the simplified tree.
        normalized = self.cancelCoefficientArtifacts(expr_id)
        if normalized != expr_id:
            candidates.append(normalized)

        search_roots = list(candidates)
        node_count = self.exprTreeNodeCount(expr_id, 5000)
        symbol_count = self.exprSymbolCost(expr_id)
        if node_count <= 4000 and symbol_count <= 900:
            for root in search_roots:
                factored = self.factorAddSubsetsRecursive(root)
                if factored != root:
                    candidates.append(factored)
                    candidates.append(self.cancelCoefficientArtifacts(factored))
                repeated = self.factorRepeatedAddSubexpr(root)
                if repeated is not None:
                    candidates.append(repeated)
                    candidates.append(self.factorAddSubsetsRecursive(repeated))
                    candidates.append(self.cancelCoefficientArtifacts(repeated))
        if node_count <= 450 and symbol_count <= 220:
            for root in list(candidates):
                expanded = self.expandAdditiveProducts(root)
                if expanded != root:
                    candidates.append(expanded)
                    candidates.append(self.factorAddSubsetsRecursive(expanded))
                    candidates.append(self.cancelCoefficientArtifacts(expanded))
                    repeated_expanded = self.factorRepeatedAddSubexpr(expanded)
                    if repeated_expanded is not None:
                        candidates.append(repeated_expanded)
                        candidates.append(self.factorAddSubsetsRecursive(repeated_expanded))

        seen: set[int] = set()
        for candidate in candidates:
            if candidate in seen:
                continue
            seen.add(candidate)
            text = self.renderExprOptimized(candidate)
            if self.betterCoefficientCandidate(text, best_text):
                best_expr = candidate
                best_text = text
        self._coefficient_post_cache[expr_id] = best_expr
        return best_expr

    def betterCoefficientCandidate(self, candidate: str, current: str) -> bool:
        if self.clearerCoefficientText(candidate, current):
            return True
        if self.clearerCoefficientText(current, candidate):
            return False
        return self.betterText(candidate, current)

    def clearerCoefficientText(self, candidate: str, current: str) -> bool:
        candidate_symbols = self.renderedSymbolCount(candidate)
        current_symbols = self.renderedSymbolCount(current)
        if candidate_symbols > current_symbols + 4:
            return False
        if self.fractionNestingPenalty(candidate) > self.fractionNestingPenalty(current):
            return False
        if len(candidate) + 20 < len(current):
            return True
        return False

    def cancelCoefficientArtifacts(self, expr_id: int) -> int:
        current = expr_id
        for _ in range(6):
            next_expr = self.absorbRationalIntoAddFactors(current)
            next_expr = self.distributeCancellableReciprocalProducts(next_expr)
            next_expr = self.factorAddSubsetsRecursive(next_expr)
            if next_expr == current:
                break
            current = next_expr
        return current

    def absorbRationalIntoAddFactors(self, expr_id: int) -> int:
        cached = self._rational_absorb_cache.get(expr_id)
        if cached is not None:
            return cached
        kind = self.expr.kind[expr_id]
        if kind in {ZERO, ONE, RATIONAL, ATOM}:
            result = expr_id
        elif kind == ADD:
            result = self.expr.add(*(self.absorbRationalIntoAddFactors(child) for child in self.expr.nodeChildren(expr_id)))
        elif kind == POW_INT:
            base, exp = self.expr.value[expr_id]
            result = self.expr.powInt(self.absorbRationalIntoAddFactors(int(base)), int(exp))
        elif kind == MUL:
            children = [self.absorbRationalIntoAddFactors(child) for child in self.expr.nodeChildren(expr_id)]
            rational_coeff = Fraction(1)
            non_rational: list[int] = []
            for child in children:
                if self.expr.kind[child] == RATIONAL:
                    p, q = self.expr.value[child]
                    rational_coeff *= Fraction(int(p), int(q))
                else:
                    non_rational.append(child)
            best = self.expr.mul(self.expr.rational(rational_coeff), *non_rational)
            if rational_coeff != 1:
                for idx, child in enumerate(non_rational):
                    if self.expr.kind[child] != ADD:
                        continue
                    common = self.commonRationalDivisor(child)
                    if common is None:
                        continue
                    p, q = self.expr.value[common]
                    common_fraction = Fraction(int(p), int(q))
                    quotient = self.expr.divByExpr(child, common)
                    if quotient is None or quotient == child:
                        continue
                    trial_factors = list(non_rational)
                    trial_factors[idx] = quotient
                    trial = self.expr.mul(self.expr.rational(rational_coeff * common_fraction), *trial_factors)
                    if self.betterExpr(trial, best) or self.betterText(self.renderExprOptimized(trial), self.renderExprOptimized(best)):
                        best = trial
            result = best
        else:
            result = expr_id
        self._rational_absorb_cache[expr_id] = result
        return result

    def distributeCancellableReciprocalProducts(self, expr_id: int) -> int:
        cached = self._cancellable_reciprocal_cache.get(expr_id)
        if cached is not None:
            return cached
        kind = self.expr.kind[expr_id]
        if kind in {ZERO, ONE, RATIONAL, ATOM}:
            result = expr_id
        elif kind == ADD:
            result = self.expr.add(*(self.distributeCancellableReciprocalProducts(child) for child in self.expr.nodeChildren(expr_id)))
        elif kind == POW_INT:
            base, exp = self.expr.value[expr_id]
            result = self.expr.powInt(self.distributeCancellableReciprocalProducts(int(base)), int(exp))
        elif kind == MUL:
            children = [self.distributeCancellableReciprocalProducts(child) for child in self.expr.nodeChildren(expr_id)]
            best = self.expr.mul(*children)
            for rec_idx, child in enumerate(children):
                if not self.isReciprocalAtomExpr(child):
                    continue
                base_expr = self.baseExprForReciprocalAtom(child)
                for add_idx, add_child in enumerate(children):
                    if add_idx == rec_idx or self.expr.kind[add_child] != ADD:
                        continue
                    terms = list(self.expr.nodeChildren(add_child))
                    if len(terms) > 96:
                        continue
                    any_cancelled = False
                    rebuilt_terms: list[int] = []
                    rest = [children[i] for i in range(len(children)) if i not in {rec_idx, add_idx}]
                    for term in terms:
                        quotient = self.expr.divByExpr(term, base_expr)
                        if quotient is None:
                            rebuilt_terms.append(self.expr.mul(*(rest + [child, term])))
                        else:
                            any_cancelled = True
                            rebuilt_terms.append(self.expr.mul(*(rest + [quotient])))
                    if not any_cancelled:
                        continue
                    candidate = self.expr.add(*rebuilt_terms)
                    candidate = self.absorbRationalIntoAddFactors(candidate)
                    candidate = self.factorAddSubsetsRecursive(candidate)
                    if self.betterExpr(candidate, best) or self.betterText(self.renderExprOptimized(candidate), self.renderExprOptimized(best)):
                        best = candidate
            result = best
        else:
            result = expr_id
        self._cancellable_reciprocal_cache[expr_id] = result
        return result

    def isReciprocalAtomExpr(self, expr_id: int) -> bool:
        if self.expr.kind[expr_id] != ATOM:
            return False
        atom_id = int(self.expr.value[expr_id])
        return self.expr.atoms.kind[atom_id] == "reciprocal_atom"

    def expandAdditiveProducts(self, expr_id: int) -> int:
        cached = self._coefficient_expand_cache.get(expr_id)
        if cached is not None:
            return cached
        if self.exprTreeNodeCount(expr_id, 5000) > 5000:
            self._coefficient_expand_cache[expr_id] = expr_id
            return expr_id
        result = self.expandAdditiveProductsInner(expr_id, 160)
        self._coefficient_expand_cache[expr_id] = result
        return result

    def expandAdditiveProductsInner(self, expr_id: int, term_limit: int) -> int:
        if self.isAggregateLinearExpr(expr_id):
            return expr_id
        kind = self.expr.kind[expr_id]
        if kind in {ZERO, ONE, RATIONAL, ATOM}:
            return expr_id
        if kind == POW_INT:
            base, exp = self.expr.value[expr_id]
            new_base = self.expandAdditiveProductsInner(int(base), term_limit)
            return self.expr.powInt(new_base, int(exp))
        if kind == ADD:
            return self.expr.add(*(self.expandAdditiveProductsInner(child, term_limit) for child in self.expr.nodeChildren(expr_id)))
        if kind != MUL:
            return expr_id
        factors = [self.expandAdditiveProductsInner(child, term_limit) for child in self.expr.nodeChildren(expr_id)]
        choices: list[tuple[int, ...]] = []
        product_terms = 1
        additive_count = 0
        for factor in factors:
            if self.expr.kind[factor] == ADD and not self.isAggregateLinearExpr(factor):
                terms = self.expr.nodeChildren(factor)
                product_terms *= len(terms)
                additive_count += 1
                if product_terms > term_limit:
                    return self.expr.mul(*factors)
                choices.append(terms)
            else:
                choices.append((factor,))
        if additive_count == 0:
            return self.expr.mul(*factors)
        terms = [self.expr.one]
        for option_group in choices:
            new_terms: list[int] = []
            for partial in terms:
                for option in option_group:
                    new_terms.append(self.expr.mul(partial, option))
            terms = new_terms
            if len(terms) > term_limit:
                return self.expr.mul(*factors)
        expanded = self.expr.add(*terms)
        if self.exprTreeNodeCount(expanded, 12000) > 12000:
            return self.expr.mul(*factors)
        return expanded

    def normalizeCoefficientExpr(self, expr_id: int) -> int:
        cached = self._coefficient_normalize_cache.get(expr_id)
        if cached is not None:
            return cached
        kind = self.expr.kind[expr_id]
        if kind in {ZERO, ONE, RATIONAL, ATOM}:
            result = expr_id
        elif kind == ADD:
            normalized_children = [self.normalizeCoefficientExpr(child) for child in self.expr.nodeChildren(expr_id)]
            result = self.expr.add(*normalized_children)
            result = self.factorCoefficientExpr(result)
        elif kind == MUL:
            normalized_children = [self.normalizeCoefficientExpr(child) for child in self.expr.nodeChildren(expr_id)]
            result = self.expr.mul(*normalized_children)
        elif kind == POW_INT:
            base, exp = self.expr.value[expr_id]
            result = self.expr.powInt(self.normalizeCoefficientExpr(int(base)), int(exp))
        else:
            result = expr_id
        self._coefficient_normalize_cache[expr_id] = result
        return result

    def factorCoefficientExpr(self, expr_id: int, allow_expand: bool = True) -> int:
        if allow_expand:
            cached = self._coefficient_factor_cache.get(expr_id)
            if cached is not None:
                return cached
        if self.expr.kind[expr_id] != ADD:
            if allow_expand:
                self._coefficient_factor_cache[expr_id] = expr_id
            return expr_id
        current = expr_id
        for _ in range(10):
            current = self.factorAddSubsets(current)
            repeated_expr = self.factorRepeatedAddSubexpr(current)
            if repeated_expr is not None and self.betterExpr(repeated_expr, current):
                current = repeated_expr
                continue
            break
        if allow_expand:
            self._coefficient_factor_cache[expr_id] = current
        return current

    def factorAddSubsets(self, expr_id: int) -> int:
        if self.expr.kind[expr_id] != ADD:
            return expr_id
        current = expr_id
        for _ in range(8):
            children = list(self.expr.nodeChildren(current))
            if len(children) < 2:
                break
            best_expr = current
            candidates = self.addSubsetDivisorCandidates(children)
            candidate_limit = 128 if len(children) <= 128 else 72
            groups = self.chooseDivisorGroups(children, candidates[:candidate_limit])
            if groups:
                consumed: set[int] = set()
                rebuilt: list[int] = []
                for divisor, indices in groups:
                    index_set = set(indices)
                    if consumed & index_set:
                        continue
                    quotients: list[int] = []
                    ok = True
                    for idx in indices:
                        quotient = self.expr.divByExpr(children[idx], divisor)
                        if quotient is None:
                            ok = False
                            break
                        quotients.append(quotient)
                    if not ok or len(quotients) < 2:
                        continue
                    consumed.update(index_set)
                    rebuilt.append(self.expr.mul(divisor, self.expr.add(*quotients)))
                rebuilt.extend(children[idx] for idx in range(len(children)) if idx not in consumed)
                candidate_expr = self.expr.add(*rebuilt)
                if candidate_expr != current and self.betterExpr(candidate_expr, best_expr):
                    best_expr = candidate_expr

            rational_divisor = self.commonRationalDivisor(current)
            if rational_divisor is not None:
                quotient = self.expr.divByExpr(current, rational_divisor)
                if quotient is not None and quotient != current:
                    candidate_expr = self.expr.mul(rational_divisor, quotient)
                    if self.betterExpr(candidate_expr, best_expr):
                        best_expr = candidate_expr

            if best_expr == current:
                break
            current = best_expr
        return current

    def betterExpr(self, candidate: int, current: int) -> bool:
        return self.exprMetric(candidate) < self.exprMetric(current)

    def exprMetric(self, expr_id: int) -> tuple[int, int, int]:
        cached = self._expr_metric_cache.get(expr_id)
        if cached is not None:
            return cached
        metric = (self.exprSymbolCost(expr_id), self.rationalSpreadCost(expr_id), self.exprTreeNodeCount(expr_id, 100000))
        self._expr_metric_cache[expr_id] = metric
        return metric

    def rationalSpreadCost(self, expr_id: int) -> int:
        kind = self.expr.kind[expr_id]
        if kind == RATIONAL:
            p, q = self.expr.value[expr_id]
            return 0 if int(q) == 1 else 1
        if kind in {ZERO, ONE, ATOM}:
            return 0
        if kind == POW_INT:
            base, exp = self.expr.value[expr_id]
            return int(exp) * self.rationalSpreadCost(int(base))
        return sum(self.rationalSpreadCost(child) for child in self.expr.nodeChildren(expr_id))

    def expandCoefficientExpr(self, expr_id: int, term_limit: int) -> int:
        key = (expr_id, term_limit)
        cached = self._expand_cache.get(key)
        if cached is not None:
            return cached
        if self.aggregateExactName(expr_id) is not None:
            self._expand_cache[key] = expr_id
            return expr_id
        kind = self.expr.kind[expr_id]
        if kind in {ZERO, ONE, RATIONAL, ATOM}:
            out = expr_id
        elif kind == ADD:
            out = self.expr.add(*(self.expandCoefficientExpr(child, term_limit) for child in self.expr.nodeChildren(expr_id)))
        elif kind == MUL:
            terms = [self.expr.one]
            for child in self.expr.nodeChildren(expr_id):
                expanded_child = self.expandCoefficientExpr(child, term_limit)
                if self.expr.kind[expanded_child] == ADD and self.aggregateExactName(expanded_child) is None:
                    child_terms = list(self.expr.nodeChildren(expanded_child))
                else:
                    child_terms = [expanded_child]
                if len(terms) * len(child_terms) > term_limit:
                    terms = [self.expr.mul(term, expanded_child) for term in terms]
                    continue
                terms = [self.expr.mul(term, child_term) for term in terms for child_term in child_terms]
            out = self.expr.add(*terms)
        elif kind == POW_INT:
            base, exp = self.expr.value[expr_id]
            expanded_base = self.expandCoefficientExpr(int(base), term_limit)
            out = self.expr.powInt(expanded_base, int(exp))
        else:
            out = expr_id
        self._expand_cache[key] = out
        return out


    def factorAddSubsetsRecursive(self, expr_id: int, depth: int = 0) -> int:
        if depth > 10:
            return expr_id
        if self.isAggregateLinearExpr(expr_id):
            return expr_id
        kind = self.expr.kind[expr_id]
        if kind in {ZERO, ONE, RATIONAL, ATOM}:
            return expr_id
        if kind == ADD:
            rebuilt = self.expr.add(*(self.factorAddSubsetsRecursive(child, depth + 1) for child in self.expr.nodeChildren(expr_id)))
            return self.factorAddSubsets(rebuilt)
        if kind == MUL:
            return self.expr.mul(*(self.factorAddSubsetsRecursive(child, depth + 1) for child in self.expr.nodeChildren(expr_id)))
        if kind == POW_INT:
            base, exp = self.expr.value[expr_id]
            return self.expr.powInt(self.factorAddSubsetsRecursive(int(base), depth + 1), int(exp))
        return expr_id

    def factorRepeatedAddSubexpr(self, expr_id: int) -> int | None:
        candidates = self.repeatedAddSubexprCandidates(expr_id)
        best_expr: int | None = None
        best_text = self.renderExprOptimized(expr_id)
        best_nodes = self.exprTreeNodeCount(expr_id, 100000)
        if self.exprTreeNodeCount(expr_id, 500) > 500 or self.exprSymbolCost(expr_id) > 120:
            return None
        for factor in candidates[:6]:
            coefficient, remainder = self.collectLinearFactor(expr_id, factor)
            if coefficient == self.expr.zero:
                continue
            if self.containsExpr(coefficient, factor):
                continue
            candidate = self.expr.add(self.expr.mul(factor, coefficient), remainder)
            if candidate == expr_id:
                continue
            candidate = self.factorAddSubsetsRecursive(candidate)
            candidate_text = self.renderExprOptimized(candidate)
            nodes = self.exprTreeNodeCount(candidate, 100000)
            if self.betterText(candidate_text, best_text) or (candidate_text == best_text and nodes < best_nodes):
                best_expr = candidate
                best_text = candidate_text
                best_nodes = nodes
        return best_expr

    def repeatedAddSubexprCandidates(self, expr_id: int) -> list[int]:
        cached = self._repeated_add_cache.get(expr_id)
        if cached is not None:
            return cached
        counts: dict[int, int] = {}

        def countOnly(node_id: int) -> None:
            if self.expr.kind[node_id] == ADD:
                cost = self.exprSymbolCost(node_id)
                if 2 <= cost <= 24:
                    counts[node_id] = counts.get(node_id, 0) + 1
            for child_id in self.expr.nodeChildren(node_id):
                countOnly(child_id)

        countOnly(expr_id)
        out = [node_id for node_id, count in counts.items() if count >= 2 and node_id != expr_id]
        out.sort(key=lambda node_id: (-(counts[node_id] - 1) * self.exprSymbolCost(node_id), -counts[node_id], -self.exprSymbolCost(node_id), self.expr.sortKey(node_id)))
        out = out[:8]
        self._repeated_add_cache[expr_id] = out
        return out

    def collectLinearFactor(self, expr_id: int, factor: int) -> tuple[int, int]:
        key = (expr_id, factor)
        cached = self._collect_linear_cache.get(key)
        if cached is not None:
            return cached
        if expr_id == self.expr.zero:
            result = (self.expr.zero, self.expr.zero)
        elif expr_id == factor:
            result = (self.expr.one, self.expr.zero)
        else:
            kind = self.expr.kind[expr_id]
            if kind in {ZERO, ONE, RATIONAL, ATOM}:
                result = (self.expr.zero, expr_id)
            elif kind == ADD:
                coefficients: list[int] = []
                remainders: list[int] = []
                for child in self.expr.nodeChildren(expr_id):
                    coefficient, remainder = self.collectLinearFactor(child, factor)
                    if coefficient != self.expr.zero:
                        coefficients.append(coefficient)
                    if remainder != self.expr.zero:
                        remainders.append(remainder)
                result = (self.expr.add(*coefficients), self.expr.add(*remainders))
            elif kind == MUL:
                children = list(self.expr.nodeChildren(expr_id))
                hit_index = -1
                for idx, child in enumerate(children):
                    if self.containsExpr(child, factor):
                        if hit_index >= 0:
                            result = (self.expr.zero, expr_id)
                            self._collect_linear_cache[key] = result
                            return result
                        hit_index = idx
                if hit_index < 0:
                    result = (self.expr.zero, expr_id)
                else:
                    child_coefficient, child_remainder = self.collectLinearFactor(children[hit_index], factor)
                    if child_coefficient == self.expr.zero:
                        result = (self.expr.zero, expr_id)
                    else:
                        rest = self.expr.mul(*(children[:hit_index] + children[hit_index + 1 :]))
                        coefficient = self.expr.mul(rest, child_coefficient)
                        remainder = self.expr.mul(rest, child_remainder) if child_remainder != self.expr.zero else self.expr.zero
                        result = (coefficient, remainder)
            else:
                result = (self.expr.zero, expr_id)
        self._collect_linear_cache[key] = result
        return result

    def containsExpr(self, expr_id: int, target: int) -> bool:
        key = (expr_id, target)
        cached = self._contains_expr_cache.get(key)
        if cached is not None:
            return cached
        result = expr_id == target or any(self.containsExpr(child, target) for child in self.expr.nodeChildren(expr_id))
        self._contains_expr_cache[key] = result
        return result

    def chooseDivisorGroups(
        self,
        children: list[int],
        candidates: list[tuple[int, tuple[int, ...]]],
    ) -> list[tuple[int, tuple[int, ...]]]:
        scored: list[tuple[int, int, int, tuple[int, ...]]] = []
        for divisor, indices in candidates:
            unique_indices = tuple(sorted(set(indices)))
            if len(unique_indices) < 2:
                continue
            divisor_cost = self.exprSymbolCost(divisor)
            saving = (len(unique_indices) - 1) * divisor_cost
            if saving <= 0 and self.expr.kind[divisor] != RATIONAL:
                continue
            # Numeric factors have zero circuit-symbol cost but still reduce
            # TeX size and tree size; keep them as secondary candidates.
            if saving <= 0:
                saving = len(unique_indices) - 1
            scored.append((saving, divisor_cost, divisor, unique_indices))
        if not scored:
            return []
        scored.sort(key=lambda item: (-item[0], -item[1], self.expr.sortKey(item[2]), item[3]))
        if len(children) <= 64 and len(scored) <= 32:
            masks = [sum(1 << idx for idx in item[3]) for item in scored]
            from functools import lru_cache

            @lru_cache(maxsize=None)
            def bestFrom(index: int, used_mask: int) -> tuple[int, int, tuple[int, ...]]:
                if index >= len(scored):
                    return (0, 0, ())
                skip_score, skip_factor_score, skip_indices = bestFrom(index + 1, used_mask)
                mask = masks[index]
                if mask & used_mask:
                    return (skip_score, skip_factor_score, skip_indices)
                take_score_tail, take_factor_tail, take_indices_tail = bestFrom(index + 1, used_mask | mask)
                take_score = scored[index][0] + take_score_tail
                take_factor_score = scored[index][1] + take_factor_tail
                take_indices = (index,) + take_indices_tail
                if take_score > skip_score:
                    return (take_score, take_factor_score, take_indices)
                if take_score == skip_score and take_factor_score > skip_factor_score:
                    return (take_score, take_factor_score, take_indices)
                if take_score == skip_score and take_factor_score == skip_factor_score and take_indices < skip_indices:
                    return (take_score, take_factor_score, take_indices)
                return (skip_score, skip_factor_score, skip_indices)

            score, _, selected = bestFrom(0, 0)
            if score <= 0:
                return []
            return [(scored[idx][2], scored[idx][3]) for idx in selected]

        selected: list[tuple[int, tuple[int, ...]]] = []
        used: set[int] = set()
        for _, _, divisor, indices in scored:
            if any(idx in used for idx in indices):
                continue
            selected.append((divisor, indices))
            used.update(indices)
        return selected

    def addSubsetDivisorCandidates(self, children: list[int]) -> list[tuple[int, tuple[int, ...]]]:
        term_items = [self.expr.factorPowerItems(child) for child in children]
        term_maps = [dict(items) for items in term_items]
        factor_to_indices: dict[int, list[int]] = {}
        for idx, items in enumerate(term_items):
            for factor, power in items:
                if power <= 0 or factor in (self.expr.zero, self.expr.one):
                    continue
                cost = self.exprSymbolCost(factor)
                if cost <= 0:
                    continue
                factor_to_indices.setdefault(factor, []).append(idx)

        out: list[tuple[int, int, tuple[int, ...]]] = []
        seen: set[tuple[int, tuple[int, ...]]] = set()

        def addCandidate(divisor: int, indices_iter) -> None:
            indices = tuple(sorted(set(indices_iter)))
            if len(indices) < 2:
                return
            if divisor in (self.expr.zero, self.expr.one):
                return
            key = (divisor, indices)
            if key in seen:
                return
            seen.add(key)
            saving = (len(indices) - 1) * self.exprSymbolCost(divisor)
            if saving <= 0:
                return
            out.append((saving, divisor, indices))

        for factor, indices_list in factor_to_indices.items():
            addCandidate(factor, indices_list)
            coeffs: list[Fraction] = []
            for idx in indices_list:
                coeff, _ = self.expr.splitCoeff(children[idx])
                if coeff != 0:
                    coeffs.append(abs(coeff))
            rational_common = self.commonRationalFromFractions(coeffs)
            if rational_common is not None and rational_common != 1:
                addCandidate(self.expr.mul(self.expr.rational(rational_common), factor), indices_list)

        rational_to_indices: dict[int, list[int]] = {}
        for idx, child in enumerate(children):
            coeff, _ = self.expr.splitCoeff(child)
            abs_coeff = abs(coeff)
            if abs_coeff == 0 or abs_coeff == 1:
                continue
            rational_to_indices.setdefault(self.expr.rational(abs_coeff), []).append(idx)
            if abs_coeff.denominator != 1:
                rational_to_indices.setdefault(self.expr.rational(Fraction(1, abs_coeff.denominator)), []).append(idx)
            if abs_coeff.numerator != 1:
                rational_to_indices.setdefault(self.expr.rational(abs_coeff.numerator), []).append(idx)
        for divisor, indices_list in rational_to_indices.items():
            if len(indices_list) >= 2:
                addCandidate(divisor, indices_list)

        # Also test common products of factors. A single-factor pass can miss
        # cases where factoring a*b is better than choosing a or b first.
        if 3 <= len(children) <= 72:
            pair_limit = min(len(children), 48)
            for left in range(pair_limit):
                left_map = term_maps[left]
                if not left_map:
                    continue
                for right in range(left + 1, pair_limit):
                    right_map = term_maps[right]
                    if not right_map:
                        continue
                    common_parts: list[int] = []
                    for factor in set(left_map).intersection(right_map):
                        power = min(left_map[factor], right_map[factor])
                        if power <= 0 or self.exprSymbolCost(factor) <= 0:
                            continue
                        common_parts.append(factor if power == 1 else self.expr.powInt(factor, power))
                    if len(common_parts) < 2:
                        continue
                    divisor = self.expr.mul(*common_parts)
                    if self.exprSymbolCost(divisor) < 2:
                        continue
                    indices = []
                    for idx, factor_map in enumerate(term_maps):
                        ok = True
                        for part in common_parts:
                            if self.expr.kind[part] == POW_INT:
                                base, exp = self.expr.value[part]
                                if factor_map.get(int(base), 0) < int(exp):
                                    ok = False
                                    break
                            elif factor_map.get(part, 0) < 1:
                                ok = False
                                break
                        if ok:
                            indices.append(idx)
                    addCandidate(divisor, indices)

        out.sort(key=lambda item: (-item[0], -len(item[2]), -self.exprSymbolCost(item[1]), self.expr.sortKey(item[1]), item[2]))
        return [(divisor, indices) for _, divisor, indices in out]

    def commonRationalFromFractions(self, coeffs: list[Fraction]) -> Fraction | None:
        if not coeffs:
            return None
        num_gcd = abs(coeffs[0].numerator)
        den_lcm = coeffs[0].denominator
        for coeff in coeffs[1:]:
            num_gcd = math.gcd(num_gcd, abs(coeff.numerator))
            den_lcm = den_lcm * coeff.denominator // math.gcd(den_lcm, coeff.denominator)
        if num_gcd == 0:
            return None
        return Fraction(num_gcd, den_lcm)

    def exprSymbolCost(self, expr_id: int) -> int:
        cached = self._expr_symbol_cost_cache.get(expr_id)
        if cached is not None:
            return cached
        aggregate_name = self.aggregate_expr_to_name.get(expr_id)
        if aggregate_name is not None and aggregate_name not in self.aggregate_disabled:
            self._expr_symbol_cost_cache[expr_id] = 1
            return 1
        kind = self.expr.kind[expr_id]
        if kind == ATOM:
            cost = 1
        elif kind in {ZERO, ONE, RATIONAL}:
            cost = 0
        elif kind == POW_INT:
            base, exp = self.expr.value[expr_id]
            cost = int(exp) * self.exprSymbolCost(int(base))
        else:
            cost = sum(self.exprSymbolCost(child) for child in self.expr.nodeChildren(expr_id))
        self._expr_symbol_cost_cache[expr_id] = cost
        return cost

    def renderExprOptimized(self, expr_id: int) -> str:
        aggregate_name = self.aggregateExactName(expr_id)
        if aggregate_name is not None:
            return aggregate_name
        cached = self._optimized_expr_cache.get(expr_id)
        if cached is not None:
            return cached
        node_count = self.exprTreeNodeCount(expr_id, 400)
        large_expr = node_count > 400
        self._optimized_expr_cache[expr_id] = self.renderExpr(expr_id)
        raw = self.renderExpr(expr_id) if large_expr else self.renderExprWithOptimizedChildren(expr_id)
        best = raw

        if not large_expr:
            factored = self.renderCommonFactorCandidate(expr_id)
            if factored is not None and self.betterText(factored, best):
                best = factored
            combined = self.renderCommonRationalReciprocalCandidate(expr_id)
            if combined is not None and self.betterText(combined, best):
                best = combined

        specs = self.reciprocalClearSpecs(expr_id)
        if large_expr:
            specs = specs[:4]
        for spec in specs:
            numerator = expr_id
            denominator = self.expr.one
            for rec_atom, count in spec:
                base_expr = self.baseExprForReciprocalAtom(rec_atom)
                for _ in range(count):
                    numerator = self.clearOneReciprocal(numerator, rec_atom)
                    denominator = self.expr.mul(denominator, base_expr)
            if large_expr:
                num_text = self.renderExpr(numerator)
                den_text = self.renderExpr(denominator)
            else:
                num_text = self.renderExprWithOptimizedChildren(numerator)
                den_text = self.renderExprWithOptimizedChildren(denominator)
            candidate = self.fractionText(num_text, den_text)
            if self.betterText(candidate, best):
                best = candidate

            if not large_expr:
                nested_factored = self.renderCommonFactorCandidate(numerator)
                if nested_factored is not None:
                    candidate = self.fractionText(nested_factored, den_text)
                    if self.betterText(candidate, best):
                        best = candidate

        self._optimized_expr_cache[expr_id] = best
        return best

    def renderCommonRationalReciprocalCandidate(self, expr_id: int) -> str | None:
        if self.expr.kind[expr_id] != ADD:
            return None
        rational_divisor = self.commonRationalDivisor(expr_id)
        if rational_divisor is None:
            return None
        quotient = self.expr.divByExpr(expr_id, rational_divisor)
        if quotient is None or quotient == expr_id:
            return None
        specs = self.reciprocalClearSpecs(quotient)
        if not specs:
            return None
        p, q = self.expr.value[rational_divisor]
        numerator_factor = self.expr.rational(int(p))
        denominator_factor = self.expr.rational(int(q))
        best: str | None = None
        for spec in specs[:8]:
            numerator = quotient
            denominator = self.expr.one
            for rec_atom, count in spec:
                base_expr = self.baseExprForReciprocalAtom(rec_atom)
                for _ in range(count):
                    numerator = self.clearOneReciprocal(numerator, rec_atom)
                    denominator = self.expr.mul(denominator, base_expr)
            numerator = self.expr.mul(numerator_factor, numerator)
            denominator = self.expr.mul(denominator_factor, denominator)
            if denominator == self.expr.one:
                continue
            numerator = self.factorAddSubsetsRecursive(numerator)
            num_text = self.renderExprWithOptimizedChildren(numerator)
            den_text = self.renderExprOptimized(denominator)
            candidate = self.fractionText(num_text, den_text)
            if best is None or self.betterText(candidate, best):
                best = candidate
        return best

    def commonRationalDivisor(self, expr_id: int) -> int | None:
        if self.expr.kind[expr_id] != ADD:
            return None
        coeffs: list[Fraction] = []
        for child in self.expr.nodeChildren(expr_id):
            coeff, _ = self.expr.splitCoeff(child)
            if coeff == 0:
                continue
            coeffs.append(abs(coeff))
        if not coeffs:
            return None
        num_gcd = abs(coeffs[0].numerator)
        den_lcm = coeffs[0].denominator
        for coeff in coeffs[1:]:
            num_gcd = math.gcd(num_gcd, abs(coeff.numerator))
            den_lcm = den_lcm * coeff.denominator // math.gcd(den_lcm, coeff.denominator)
        common = Fraction(num_gcd, den_lcm)
        if common == 1:
            return None
        return self.expr.rational(common)

    def factorQuotientText(self, divisor: int, quotient: int) -> str:
        quotient_text = self.renderExpr(quotient) if self.exprTreeNodeCount(quotient, 120) > 120 else self.renderExprWithOptimizedChildren(quotient)
        if self.expr.kind[divisor] == RATIONAL:
            p, q = self.expr.value[divisor]
            frac = Fraction(int(p), int(q))
            sign = "-" if frac < 0 else ""
            abs_frac = abs(frac)
            if abs_frac.numerator == 1:
                best = sign + self.fractionText(quotient_text, str(abs_frac.denominator))
                forward = sign + self.productText(self.fractionText("1", str(abs_frac.denominator)), quotient_text)
                if self.longSimpleDenominatorPenalty(best) > 0 or self.betterText(forward, best):
                    best = forward
            else:
                numerator_text = self.productText(str(abs_frac.numerator), quotient_text)
                best = sign + self.fractionText(numerator_text, str(abs_frac.denominator))
                forward = sign + self.productText(self.fractionText(str(abs_frac.numerator), str(abs_frac.denominator)), quotient_text)
                if self.longSimpleDenominatorPenalty(best) > 0 or self.betterText(forward, best):
                    best = forward
            distributed = self.renderDistributedRationalProduct(frac, quotient)
            if distributed is not None and self.betterText(distributed, best):
                return distributed
            return best
        return self.productText(self.renderExprOptimized(divisor), quotient_text)

    def renderDistributedRationalProduct(self, coeff: Fraction, expr_id: int) -> str | None:
        if self.expr.kind[expr_id] != ADD:
            return None
        terms: list[str] = []
        coeff_expr = self.expr.rational(coeff)
        for child in self.expr.nodeChildren(expr_id):
            term_expr = self.expr.mul(coeff_expr, child)
            terms.append(self.renderExprOptimized(term_expr))
        return self.joinSignedTerms(terms)

    def renderCommonFactorCandidate(self, expr_id: int) -> str | None:
        if self.expr.kind[expr_id] != ADD:
            return None
        if self.exprTreeNodeCount(expr_id, 420) > 420:
            return None
        factors = self.expr.commonTermFactorPowerItems(expr_id)
        if not factors:
            return None
        candidates: list[int] = []
        rational_divisor = self.commonRationalDivisor(expr_id)
        if rational_divisor is not None:
            candidates.append(rational_divisor)
        for factor, power in factors:
            if factor in (self.expr.zero, self.expr.one) or power <= 0:
                continue
            candidates.append(self.expr.powInt(factor, power) if power > 1 else factor)
        if not candidates:
            return None
        candidates.sort(key=lambda item: (-len(self.renderExpr(item)), item))
        if len(candidates) > 1:
            candidates.insert(0, self.expr.mul(*candidates))
        best: str | None = None
        for divisor in candidates[:10]:
            quotient = self.expr.divByExpr(expr_id, divisor)
            if quotient is None or quotient == expr_id:
                continue
            candidate = self.factorQuotientText(divisor, quotient)
            if best is None or self.betterText(candidate, best):
                best = candidate
        return best

    def renderExprWithOptimizedChildren(self, expr_id: int) -> str:
        kind = self.expr.kind[expr_id]
        if kind in {ZERO, ONE, RATIONAL, ATOM}:
            return self.renderExpr(expr_id)
        if kind == ADD:
            return self.renderAddWithRenderer(expr_id, self.renderExprOptimized)
        if kind == MUL:
            return self.renderMulWithRenderer(expr_id, self.renderExprOptimized)
        if kind == POW_INT:
            base, exp = self.expr.value[expr_id]
            base_text = self.renderExprOptimized(int(base))
            if self.needsParens(base_text):
                base_text = r"\left(" + base_text + r"\right)"
            return base_text + f"^{{{int(exp)}}}"
        return self.renderExpr(expr_id)

    def reciprocalMaxPower(self, expr_id: int, atom_expr: int) -> int:
        kind = self.expr.kind[expr_id]
        if expr_id == atom_expr:
            return 1
        if kind in {ZERO, ONE, RATIONAL, ATOM}:
            return 0
        if kind == ADD:
            return max((self.reciprocalMaxPower(child, atom_expr) for child in self.expr.nodeChildren(expr_id)), default=0)
        if kind == MUL:
            return sum(self.reciprocalMaxPower(child, atom_expr) for child in self.expr.nodeChildren(expr_id))
        if kind == POW_INT:
            base, exp = self.expr.value[expr_id]
            return int(exp) * self.reciprocalMaxPower(int(base), atom_expr)
        return sum(self.reciprocalMaxPower(child, atom_expr) for child in self.expr.nodeChildren(expr_id))

    def reciprocalClearSpecs(self, expr_id: int) -> list[tuple[tuple[int, int], ...]]:
        atoms = self.reciprocalAtomExprs(expr_id)
        if not atoms:
            return []
        limited = atoms[:8]
        specs: list[tuple[tuple[int, int], ...]] = []
        seen: set[tuple[tuple[int, int], ...]] = set()

        def addSpec(items: list[tuple[int, int]]) -> None:
            spec = tuple((atom, count) for atom, count in items if count > 0)
            if spec and spec not in seen:
                seen.add(spec)
                specs.append(spec)

        # Common-denominator candidates are tried first.  The count is the
        # maximum reciprocal power in any additive term, not the total number
        # of appearances, so a/r+b/r becomes (a+b)/r rather than an over-cleared
        # expression.  This matches the relaxed final-E rational form.
        max_power_items = [(atom, self.reciprocalMaxPower(expr_id, atom)) for atom in limited]
        if len(max_power_items) <= 8:
            addSpec(max_power_items)
        if len(max_power_items) >= 2:
            addSpec(max_power_items[:4])
            for i in range(min(6, len(max_power_items))):
                for j in range(i + 1, min(6, len(max_power_items))):
                    addSpec([max_power_items[i], max_power_items[j]])

        for atom in limited:
            occ = self.reciprocalOccurrence(expr_id, atom)
            max_power = self.reciprocalMaxPower(expr_id, atom)
            counts = sorted(set([1, max_power, min(occ, 2), min(occ, 3), min(occ, 4), occ]))
            for count in counts:
                if count > 0:
                    addSpec([(atom, count)])

        prefix = limited[:4]
        addSpec([(atom, 1) for atom in prefix])
        addSpec([(atom, self.reciprocalMaxPower(expr_id, atom)) for atom in prefix])
        for i in range(len(prefix)):
            for j in range(i + 1, len(prefix)):
                a = prefix[i]
                b = prefix[j]
                addSpec([(a, self.reciprocalMaxPower(expr_id, a)), (b, self.reciprocalMaxPower(expr_id, b))])
        return specs

    def fractionText(self, numerator: str, denominator: str) -> str:
        num = numerator
        den = denominator
        negative = False
        if num.startswith("-") and not self.needsParens(num[1:]):
            negative = not negative
            num = num[1:]
        if den.startswith("-") and not self.needsParens(den[1:]):
            negative = not negative
            den = den[1:]
        text = r"\frac{" + num + "}{" + den + "}"
        return "-" + text if negative else text


    def exprTreeNodeCount(self, expr_id: int, limit: int) -> int:
        cached = self._expr_cost_cache.get(expr_id)
        if cached is not None:
            return cached
        kind = self.expr.kind[expr_id]
        if kind in {ZERO, ONE, RATIONAL, ATOM}:
            out = 1
        else:
            total = 1
            for child in self.expr.nodeChildren(expr_id):
                total += self.exprTreeNodeCount(child, limit)
                if total > limit:
                    self._expr_cost_cache[expr_id] = total
                    return total
            out = total
        self._expr_cost_cache[expr_id] = out
        return out

    def reciprocalAtomExprs(self, expr_id: int) -> list[int]:
        out: set[int] = set()
        self.collectReciprocalAtomExprs(expr_id, out)
        return sorted(out, key=lambda atom_expr: -self.reciprocalOccurrence(expr_id, atom_expr))

    def collectReciprocalAtomExprs(self, expr_id: int, out: set[int]) -> None:
        if self.isAggregateLinearExpr(expr_id):
            return
        kind = self.expr.kind[expr_id]
        if kind == ATOM:
            atom_id = int(self.expr.value[expr_id])
            if self.expr.atoms.kind[atom_id] == "reciprocal_atom":
                out.add(expr_id)
            return
        if kind in {ZERO, ONE, RATIONAL}:
            return
        for child in self.expr.nodeChildren(expr_id):
            self.collectReciprocalAtomExprs(child, out)

    def reciprocalOccurrence(self, expr_id: int, atom_expr: int) -> int:
        if self.isAggregateLinearExpr(expr_id):
            return 0
        kind = self.expr.kind[expr_id]
        if expr_id == atom_expr:
            return 1
        if kind in {ZERO, ONE, RATIONAL, ATOM}:
            return 0
        if kind == POW_INT:
            base, exp = self.expr.value[expr_id]
            return int(exp) if int(base) == atom_expr else self.reciprocalOccurrence(int(base), atom_expr) * int(exp)
        return sum(self.reciprocalOccurrence(child, atom_expr) for child in self.expr.nodeChildren(expr_id))



    def baseExprForReciprocalAtom(self, atom_expr: int) -> int:
        atom_id = int(self.expr.value[atom_expr])
        base_id = self.expr.atoms.base_param[atom_id]
        if base_id is None:
            raise KernelTexError("reciprocal atom has no base parameter")
        return self.expr.atom(base_id)

    def containsReciprocal(self, expr_id: int, atom_expr: int) -> bool:
        if self.isAggregateLinearExpr(expr_id):
            return False
        key = (expr_id, atom_expr)
        cached = self._contains_rec_cache.get(key)
        if cached is not None:
            return cached
        if expr_id == atom_expr:
            out = True
        else:
            kind = self.expr.kind[expr_id]
            if kind in {ZERO, ONE, RATIONAL, ATOM}:
                out = False
            else:
                out = any(self.containsReciprocal(child, atom_expr) for child in self.expr.nodeChildren(expr_id))
        self._contains_rec_cache[key] = out
        return out

    def clearOneReciprocal(self, expr_id: int, atom_expr: int) -> int:
        if self.isAggregateLinearExpr(expr_id):
            return self.expr.mul(self.baseExprForReciprocalAtom(atom_expr), expr_id)
        key = (expr_id, atom_expr)
        cached = self._clear_rec_cache.get(key)
        if cached is not None:
            return cached
        base_expr = self.baseExprForReciprocalAtom(atom_expr)
        kind = self.expr.kind[expr_id]
        if expr_id == atom_expr:
            out = self.expr.one
        elif kind == ZERO:
            out = self.expr.zero
        elif kind in {ONE, RATIONAL, ATOM}:
            out = self.expr.mul(base_expr, expr_id)
        elif kind == ADD:
            out = self.expr.add(*(self.clearOneReciprocal(child, atom_expr) for child in self.expr.nodeChildren(expr_id)))
        elif kind == MUL:
            children = list(self.expr.nodeChildren(expr_id))
            replaced = False
            for idx, child in enumerate(children):
                if self.containsReciprocal(child, atom_expr):
                    children[idx] = self.clearOneReciprocal(child, atom_expr)
                    replaced = True
                    break
            if replaced:
                out = self.expr.mul(*children)
            else:
                out = self.expr.mul(base_expr, expr_id)
        elif kind == POW_INT:
            base, exp = self.expr.value[expr_id]
            if int(base) == atom_expr and int(exp) > 0:
                out = self.expr.powInt(atom_expr, int(exp) - 1)
            elif self.containsReciprocal(int(base), atom_expr):
                # R * f^n = (R*f) * f^(n-1). This keeps the transformation exact without expanding powers.
                out = self.expr.mul(self.clearOneReciprocal(int(base), atom_expr), self.expr.powInt(int(base), int(exp) - 1))
            else:
                out = self.expr.mul(base_expr, expr_id)
        else:
            out = self.expr.mul(base_expr, expr_id)
        self._clear_rec_cache[key] = out
        return out

    def renderedSymbolCount(self, text: str) -> int:
        cleaned = re.sub(r"\\frac|\\left|\\right|\\sqrt|\\text|\\begin|\\end", " ", text)
        tokens = re.findall(r"[A-Za-z]+(?:_\{[^}]+\})?", cleaned)
        ignored = {"s", "H", "E", "Z", "P"}
        return sum(1 for token in tokens if token.split("_")[0] not in ignored)

    def renderExpr(self, expr_id: int) -> str:
        aggregate_name = self.aggregateExactName(expr_id)
        if aggregate_name is not None:
            return aggregate_name
        cached = self._expr_cache.get(expr_id)
        if cached is not None:
            return cached
        kind = self.expr.kind[expr_id]
        if kind == ZERO:
            out = "0"
        elif kind == ONE:
            out = "1"
        elif kind == RATIONAL:
            p, q = self.expr.value[expr_id]
            if int(q) == 1:
                out = str(int(p))
            else:
                sign = "-" if int(p) < 0 else ""
                out = sign + r"\frac{" + str(abs(int(p))) + "}{" + str(int(q)) + "}"
        elif kind == ATOM:
            out = self.expr.atoms.atomLatex(int(self.expr.value[expr_id]))
        elif kind == ADD:
            out = self.renderAdd(expr_id)
        elif kind == MUL:
            out = self.renderMul(expr_id)
        elif kind == POW_INT:
            base, exp = self.expr.value[expr_id]
            base_text = self.renderExpr(int(base))
            if self.needsParens(base_text):
                base_text = r"\left(" + base_text + r"\right)"
            out = base_text + f"^{{{int(exp)}}}"
        else:
            raise KernelTexError("unknown expression kind")
        self._expr_cache[expr_id] = out
        return out

    def aggregateExactName(self, expr_id: int) -> str | None:
        name = self.aggregate_expr_to_name.get(expr_id)
        if name is None or name in self.aggregate_disabled:
            return None
        self._aggregate_use_counts[name] = self._aggregate_use_counts.get(name, 0) + 1
        return name

    def betterText(self, candidate: str, current: str) -> bool:
        candidate_cost = self.renderedSymbolCount(candidate)
        current_cost = self.renderedSymbolCount(current)
        if candidate_cost != current_cost:
            return candidate_cost < current_cost
        candidate_penalty = self.fractionNestingPenalty(candidate)
        current_penalty = self.fractionNestingPenalty(current)
        if candidate_penalty != current_penalty:
            return candidate_penalty < current_penalty
        return len(candidate) < len(current)

    def fractionNestingPenalty(self, text: str) -> int:
        return (
            4 * text.count(r"\frac{1}{2}")
            + 3 * len(re.findall(r"\}\{2(?:\s|\}|[A-Za-z_])", text))
            + 2 * text.count(r"\frac{\frac")
            + 2 * text.count(r"\frac{-\frac")
            + 2 * text.count(r"}{\frac")
            + 2 * text.count(r"}{-\frac")
        )

    def renderAdd(self, expr_id: int) -> str:
        return self.renderAddWithRenderer(expr_id, self.renderExpr)

    def renderAddWithRenderer(self, expr_id: int, render_child) -> str:
        grouped = self.renderAddWithAggregateGroups(expr_id, render_child)
        if grouped is not None:
            return grouped
        return self.joinSignedTerms([render_child(child) for child in self.expr.nodeChildren(expr_id)])

    def renderAddWithAggregateGroups(self, expr_id: int, render_child) -> str | None:
        if not self.aggregate_name_to_atom_exprs:
            return None
        children = list(self.expr.nodeChildren(expr_id))
        candidates: list[tuple[int, int, str, Fraction, tuple[int, ...]]] = []
        for name in sorted(self.aggregate_name_to_atom_exprs, key=lambda item: (self.aggregate_name_to_atom_count[item], item)):
            if name in self.aggregate_disabled:
                continue
            atom_set = set(self.aggregate_name_to_atom_exprs[name])
            by_coeff: dict[Fraction, list[tuple[int, int]]] = {}
            for idx, child in enumerate(children):
                coeff, base = self.expr.splitCoeff(child)
                if base in atom_set:
                    by_coeff.setdefault(coeff, []).append((idx, base))
            for coeff, matches in by_coeff.items():
                if {base for _, base in matches} != atom_set:
                    continue
                indices = tuple(sorted(idx for idx, _ in matches))
                if len(indices) != len(atom_set):
                    continue
                # Circuit-symbol metric: replacing all members of the set by one
                # aggregate saves one symbol for every removed atom occurrence.
                benefit = len(atom_set) - 1
                if benefit <= 0:
                    continue
                candidates.append((benefit, len(indices), name, coeff, indices))
        if not candidates:
            return None
        chosen = self.chooseAggregateAddGroups(candidates, len(children))
        if not chosen:
            return None
        consumed: set[int] = set()
        replacement_texts: list[str] = []
        for _, _, name, coeff, indices in chosen:
            if any(idx in consumed for idx in indices):
                continue
            consumed.update(indices)
            replacement_texts.append(self.aggregateTermText(name, coeff))
            self._aggregate_use_counts[name] = self._aggregate_use_counts.get(name, 0) + 1
        remaining = [render_child(child) for idx, child in enumerate(children) if idx not in consumed]
        raw = self.joinSignedTerms([render_child(child) for child in children])
        grouped = self.joinSignedTerms(replacement_texts + remaining)
        return grouped if self.betterText(grouped, raw) else None

    def chooseAggregateAddGroups(
        self,
        candidates: list[tuple[int, int, str, Fraction, tuple[int, ...]]],
        child_count: int,
    ) -> list[tuple[int, int, str, Fraction, tuple[int, ...]]]:
        # Exact memoized DP for the small local set-packing problem. For very
        # large additions, use the same objective with deterministic greedy
        # selection to avoid exponential behavior.
        candidates = sorted(
            candidates,
            key=lambda item: (-item[0], -item[1], item[2], item[3].numerator, item[3].denominator, item[4]),
        )
        if child_count <= 48 and len(candidates) <= 24:
            masks = [sum(1 << idx for idx in item[4]) for item in candidates]
            from functools import lru_cache

            @lru_cache(maxsize=None)
            def bestFrom(index: int, used_mask: int) -> tuple[int, tuple[int, ...]]:
                if index >= len(candidates):
                    return (0, ())
                skip_score, skip_items = bestFrom(index + 1, used_mask)
                mask = masks[index]
                if mask & used_mask:
                    return (skip_score, skip_items)
                take_tail_score, take_tail_items = bestFrom(index + 1, used_mask | mask)
                take_score = candidates[index][0] + take_tail_score
                take_items = (index,) + take_tail_items
                if take_score > skip_score:
                    return (take_score, take_items)
                if take_score == skip_score and take_items < skip_items:
                    return (take_score, take_items)
                return (skip_score, skip_items)

            score, selected_indices = bestFrom(0, 0)
            if score <= 0:
                return []
            return [candidates[idx] for idx in selected_indices]

        selected: list[tuple[int, int, str, Fraction, tuple[int, ...]]] = []
        used: set[int] = set()
        for candidate in candidates:
            indices = candidate[4]
            if any(idx in used for idx in indices):
                continue
            selected.append(candidate)
            used.update(indices)
        return selected

    def aggregateTermText(self, name: str, coeff: Fraction) -> str:
        if coeff == 1:
            return name
        if coeff == -1:
            return "-" + name
        coeff_text = self.renderExpr(self.expr.rational(coeff))
        text = self.productText(coeff_text, name)
        return text

    def joinSignedTerms(self, term_texts: list[str]) -> str:
        if not term_texts:
            return "0"
        negative_terms = [text[1:] for text in term_texts if text.startswith("-")]
        positive_terms = [text for text in term_texts if not text.startswith("-")]
        if len(negative_terms) == len(term_texts) and len(term_texts) > 1:
            inner = self.joinSignedTerms(negative_terms)
            return r"-\left(" + inner + r"\right)"
        if positive_terms and term_texts[0].startswith("-"):
            term_texts = positive_terms + ["-" + text for text in negative_terms]
        parts: list[str] = []
        for text in term_texts:
            if text.startswith("-"):
                op = "-"
                body = text[1:]
            else:
                op = "+"
                body = text
            if not parts:
                parts.append(("-" if op == "-" else "") + body)
            else:
                parts.append(op + body)
        return "".join(parts)

    def renderMul(self, expr_id: int) -> str:
        return self.renderMulWithRenderer(expr_id, self.renderExpr)

    def renderMulWithRenderer(self, expr_id: int, render_child) -> str:
        children = list(self.expr.nodeChildren(expr_id))
        sign = ""
        factors: list[str] = []
        for child in children:
            text = render_child(child)
            if text == "-1":
                sign = "-" if not sign else ""
                continue
            if text.startswith("-") and self.expr.kind[child] == RATIONAL:
                sign = "-" if not sign else ""
                text = text[1:]
            if self.needsParens(text):
                text = r"\left(" + text + r"\right)"
            factors.append(text)
        product = sign + (" ".join(factors) if factors else "1")
        fraction = self.renderMulAsSingleFraction(children, render_child)
        if fraction is not None:
            if sign and not fraction.startswith("-"):
                fraction = sign + fraction
            if self.betterText(fraction, product):
                return fraction
        return product

    def renderMulAsSingleFraction(self, children: list[int], render_child) -> str | None:
        numerator: list[int] = []
        denominator: list[int] = []
        sign = 1
        rational_num = 1
        rational_den = 1
        for child in children:
            if self.expr.kind[child] == RATIONAL:
                p, q = self.expr.value[child]
                frac = Fraction(int(p), int(q))
                if frac < 0:
                    sign *= -1
                    frac = -frac
                rational_num *= frac.numerator
                rational_den *= frac.denominator
                continue
            if self.expr.kind[child] == ATOM:
                atom_id = int(self.expr.value[child])
                if self.expr.atoms.kind[atom_id] == "reciprocal_atom":
                    denominator.append(self.baseExprForReciprocalAtom(child))
                    continue
            if self.expr.kind[child] == POW_INT:
                base, exp = self.expr.value[child]
                if self.expr.kind[int(base)] == ATOM:
                    atom_id = int(self.expr.value[int(base)])
                    if self.expr.atoms.kind[atom_id] == "reciprocal_atom":
                        denominator.append(self.expr.powInt(self.baseExprForReciprocalAtom(int(base)), int(exp)))
                        continue
            numerator.append(child)
        if rational_den != 1:
            denominator.insert(0, self.expr.rational(rational_den))
        if rational_num != 1:
            numerator.insert(0, self.expr.rational(rational_num))
        if not denominator:
            return None
        num_expr = self.expr.mul(*numerator) if numerator else self.expr.one
        if sign < 0:
            num_expr = self.expr.neg(num_expr)
        den_expr = self.expr.mul(*denominator)
        rational_divisor = self.commonRationalDivisor(num_expr) if self.expr.kind[num_expr] == ADD else None
        if rational_divisor is not None:
            quotient = self.expr.divByExpr(num_expr, rational_divisor)
            if quotient is not None and quotient != num_expr:
                p, q = self.expr.value[rational_divisor]
                num_expr = self.expr.mul(self.expr.rational(int(p)), quotient)
                den_expr = self.expr.mul(den_expr, self.expr.rational(int(q)))
        negative_fraction = False
        if self.addTermsAllNegative(num_expr):
            negative_fraction = True
            num_expr = self.negateAddTerms(num_expr)
        num_text = render_child(num_expr)
        den_text = self.renderExprOptimized(den_expr)
        text = self.fractionText(num_text, den_text)
        if negative_fraction and not text.startswith("-"):
            text = "-" + text
        return text

    def addTermsAllNegative(self, expr_id: int) -> bool:
        if self.expr.kind[expr_id] != ADD:
            return False
        return all(self.expr.splitCoeff(child)[0] < 0 for child in self.expr.nodeChildren(expr_id))

    def negateAddTerms(self, expr_id: int) -> int:
        if self.expr.kind[expr_id] != ADD:
            return self.expr.neg(expr_id)
        return self.expr.add(*(self.expr.neg(child) for child in self.expr.nodeChildren(expr_id)))


    def needsParens(self, text: str) -> bool:
        depth = 0
        idx = 0
        while idx < len(text):
            if text.startswith(r"\left", idx):
                depth += 1
                idx += 5
                continue
            if text.startswith(r"\right)", idx):
                depth = max(0, depth - 1)
                idx += 6
                continue
            char = text[idx]
            if char == "{":
                depth += 1
            elif char == "}":
                depth = max(0, depth - 1)
            elif depth == 0 and idx > 0 and char in "+-":
                return True
            idx += 1
        return False
