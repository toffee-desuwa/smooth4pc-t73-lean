#!/usr/bin/env python3
"""Recompute the T73 divided cubic from a compact, result-free input.

The input contains primitive crossing rows and geometric chronology, not the
derived B44/B88 words and not the expected cubic.  This script independently
rebuilds both words, evaluates the unreduced Burau action in
Z[epsilon]/(epsilon^7), substitutes t=q^-2 with q=1+h, and prints a
path-independent receipt.

The cup vector u and the cap row ell are not read from hand-written lists.
They are derived from the single endpoint authority file
data/T73_ENDPOINT_CONVENTION.json (physical endpoint identities, orientations,
geometric and public orders, pivotal coefficients) by
build_t73_endpoint_transport.derive_endpoint_terms; the public input records
only the name and SHA-256 of that authority file.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
DEFAULT_INPUT = HERE.parent / "data" / "T73_DELTA3_PUBLIC_INPUT.json"
DEFAULT_RECEIPT = HERE.parent / "data" / "T73_DELTA3_PUBLIC_RECEIPT.json"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest().upper()


def canonical_sha(value: Any) -> str:
    raw = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return sha256_bytes(raw)


def inverse_word(word: list[int]) -> list[int]:
    return [-letter for letter in reversed(word)]


def left_pure_generator(i: int, j: int, sign: int) -> list[int]:
    """Right label ``j`` moves left under label ``i``."""

    if not (1 <= i < j and sign in (-1, 1)):
        raise ValueError(f"invalid L generator: i={i}, j={j}, sign={sign}")
    return (
        list(range(j - 1, i, -1))
        + [sign * i, sign * i]
        + [-k for k in range(i + 1, j)]
    )


def right_pure_generator(i: int, j: int, sign: int) -> list[int]:
    """Left label ``i`` moves right over label ``j``."""

    if not (1 <= i < j and sign in (-1, 1)):
        raise ValueError(f"invalid R generator: i={i}, j={j}, sign={sign}")
    return (
        list(range(i, j - 1))
        + [sign * (j - 1), sign * (j - 1)]
        + [-k for k in range(j - 2, i - 1, -1)]
    )


def parse_crossing_rows(data: dict[str, Any]) -> list[dict[str, Any]]:
    columns = data["point_push"]["crossing_row_columns"]
    expected_columns = [
        "source_index",
        "x",
        "y",
        "sign",
        "under_owner",
        "under_segment",
        "over_owner",
        "over_segment",
        "moving_wicket",
        "other_wicket",
        "m2_passage_exponent",
        "geometry",
    ]
    if columns != expected_columns:
        raise ValueError("unexpected crossing-row schema")
    rows = []
    for raw in data["point_push"]["crossing_rows"]:
        if len(raw) != len(columns):
            raise ValueError("malformed crossing row")
        row = dict(zip(columns, raw))
        if row["sign"] not in (-1, 1) or row["geometry"] not in ("L", "R"):
            raise ValueError("invalid crossing sign or geometry")
        if row["over_owner"] != "r_xy" or row["under_owner"] != "m_2":
            raise ValueError("crossing row is outside the frozen r_xy/m_2 sector")
        if row["geometry"] == "L":
            if row["other_wicket"] != 1 or row["moving_wicket"] <= 1:
                raise ValueError("invalid L wicket roles")
        else:
            if row["moving_wicket"] != 1 or row["other_wicket"] <= 1:
                raise ValueError("invalid R wicket roles")
        rows.append(row)
    return rows


def row_artin_word(row: dict[str, Any]) -> list[int]:
    i = min(row["moving_wicket"], row["other_wicket"])
    j = max(row["moving_wicket"], row["other_wicket"])
    if row["geometry"] == "L":
        return left_pure_generator(i, j, row["sign"])
    return right_pure_generator(i, j, row["sign"])


def build_oriented_b44(data: dict[str, Any]) -> tuple[list[int], list[int]]:
    rows = parse_crossing_rows(data)
    row_digest = canonical_sha(data["point_push"]["crossing_rows"])
    if row_digest != data["point_push"]["crossing_rows_sha256"]:
        raise ValueError("crossing-row SHA mismatch")

    # Reconstruct point-push order from the public component chronology.  This
    # is intentionally not the source emitter order: each route leg declares
    # both its r_xy segment and its oriented x direction.
    chronology = data["point_push"]["chronology"]
    if chronology["emitter_sort"] != "increasing_x_within_each_horizontal_segment":
        raise ValueError("unexpected emitter ordering rule")
    ordered_rows: list[dict[str, Any]] = []
    seen_source_indices: set[int] = set()
    for leg in chronology["oriented_horizontal_legs"]:
        direction = leg["x_direction"]
        if direction not in ("increasing", "decreasing"):
            raise ValueError("invalid oriented x direction")
        selected = [
            row
            for row in rows
            if row["geometry"] == leg["geometry"]
            and row["over_segment"] == leg["over_segment"]
        ]
        selected.sort(key=lambda row: (row["x"], row["source_index"]))
        if direction == "decreasing":
            selected.reverse()
        for row in selected:
            if row["source_index"] in seen_source_indices:
                raise ValueError("chronology selects one crossing more than once")
            seen_source_indices.add(row["source_index"])
        ordered_rows.extend(selected)
    if len(ordered_rows) != len(rows):
        raise ValueError("chronology does not cover every primitive crossing row")
    source_indices = [row["source_index"] for row in ordered_rows]
    if source_indices != data["point_push"]["oriented_source_indices"]:
        raise ValueError("oriented source-index order mismatch")

    word = [letter for row in ordered_rows for letter in row_artin_word(row)]
    expected = data["point_push"]["derived_integrity"]
    if len(ordered_rows) != expected["factor_count"]:
        raise ValueError("factor-count mismatch")
    if len(word) != expected["B44_length"]:
        raise ValueError("B44 length mismatch")
    if canonical_sha(word) != expected["B44_sha256"]:
        raise ValueError("B44 SHA mismatch")
    return word, source_indices


def cable_letter(letter: int) -> list[int]:
    i = abs(letter)
    block = [2 * i, 2 * i + 1, 2 * i - 1, 2 * i]
    return block if letter > 0 else inverse_word(block)


def cable_word(word: list[int]) -> list[int]:
    return [cabled for letter in word for cabled in cable_letter(letter)]


def poly_add(left: list[int], right: list[int]) -> list[int]:
    return [x + y for x, y in zip(left, right)]


def poly_sub(left: list[int], right: list[int]) -> list[int]:
    return [x - y for x, y in zip(left, right)]


def mul_epsilon(poly: list[int]) -> list[int]:
    return [0] + poly[:-1]


def mul_t_inverse(poly: list[int]) -> list[int]:
    # (1+epsilon)^-1 = 1-epsilon+epsilon^2-... .
    return [
        sum(((-1) ** r) * poly[k - r] for r in range(k + 1))
        for k in range(len(poly))
    ]


def apply_generator_left(vector: list[list[int]], letter: int) -> None:
    i = abs(letter) - 1
    if i < 0 or i + 1 >= len(vector):
        raise ValueError(f"Artin generator out of range: {letter}")
    x = vector[i][:]
    y = vector[i + 1][:]
    if letter > 0:
        vector[i] = poly_add(y, poly_sub(mul_epsilon(y), mul_epsilon(x)))
        vector[i + 1] = x
    else:
        t_inverse_y = mul_t_inverse(y)
        vector[i] = y
        vector[i + 1] = poly_add(mul_t_inverse(x), mul_epsilon(t_inverse_y))


def apply_word(word: list[int], vector: list[list[int]]) -> list[list[int]]:
    result = [poly[:] for poly in vector]
    # rho(letter_1 ... letter_m) acts on a column from right to left.
    for letter in reversed(word):
        apply_generator_left(result, letter)
    return result


def delta_apply(word: list[int], vector: list[list[int]]) -> list[list[int]]:
    acted = apply_word(word, vector)
    return [poly_sub(after, before) for after, before in zip(acted, vector)]


def sparse_vector(
    dimension: int, degree: int, terms: list[list[int]]
) -> list[list[int]]:
    zero = [0] * (degree + 1)
    vector = [zero[:] for _ in range(dimension)]
    for index, coefficient in terms:
        if not 0 <= index < dimension:
            raise ValueError(f"vector index out of range: {index}")
        vector[index][0] += coefficient
    return vector


def apply_covector(vector: list[list[int]], terms: list[list[int]]) -> list[int]:
    result = [0] * len(vector[0])
    for index, coefficient in terms:
        result = poly_add(result, [coefficient * value for value in vector[index]])
    return result


def poly_mul(left: list[int], right: list[int], degree: int) -> list[int]:
    result = [0] * (degree + 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            if i + j <= degree:
                result[i + j] += x * y
    return result


def epsilon_of_h(degree: int) -> list[int]:
    # epsilon=(1+h)^-2-1; coefficient of h^k is (-1)^k*(k+1).
    return [0] + [((-1) ** k) * (k + 1) for k in range(1, degree + 1)]


def substitute_epsilon_with_h(coeffs: list[int], degree: int) -> list[int]:
    epsilon = epsilon_of_h(degree)
    result = [0] * (degree + 1)
    power = [1] + [0] * degree
    for coefficient in coeffs:
        result = poly_add(result, [coefficient * value for value in power])
        power = poly_mul(power, epsilon, degree)
    return result


def free_reduce(tokens: list[str]) -> list[str]:
    inverse = {"W": "w", "w": "W", "F": "f", "f": "F"}
    stack: list[str] = []
    for token in tokens:
        if stack and inverse.get(token) == stack[-1]:
            stack.pop()
        else:
            stack.append(token)
    return stack


def check_binding(data: dict[str, Any]) -> None:
    binding = data["hattori_binding"]
    source_extract = binding["source_extract"]
    if canonical_sha(source_extract) != binding["source_extract_sha256"]:
        raise ValueError("Hattori source-extract SHA mismatch")
    required_components = {
        "eta_Rw[T0] -> Id_(W U1) tensor X^227",
        "eta_Rw[T1] -> Id_(W W^-1 U1)=Id_U1 tensor X^227",
    }
    if not required_components.issubset(source_extract["Hattori_components"]):
        raise ValueError("Hattori source extract omits a required component")
    required_objects = {
        "T0=F_Omega^-1 U1",
        "T1=F_Omega^-1 W^-1 U1",
    }
    if not required_objects.issubset(
        source_extract["framing_corrected_objects_for_B_equals_W_FOmega"]
    ):
        raise ValueError("Hattori source extract omits a required object")
    if source_extract["same_fixed_coefficient"] != "R_w[Omega]":
        raise ValueError("unexpected fixed Hattori coefficient")
    if free_reduce(binding["coefficient_word"] + binding["T0_word"]) != ["W", "U"]:
        raise ValueError("B*T0 order mismatch")
    if free_reduce(binding["coefficient_word"] + binding["T1_word"]) != ["U"]:
        raise ValueError("B*T1 order mismatch")
    if not binding["same_fixed_coefficient"]:
        raise ValueError("Hattori inputs do not use one fixed coefficient")


def derive_endpoint_terms_from_authority(input_path: Path, data: dict[str, Any]) -> dict[str, Any]:
    """Cup/cap constant terms from the single endpoint authority file.

    The terms are recomputed from the recorded transport monomials and the
    oriented coevaluation / pivotal evaluation; a hand-written ``u_terms`` or
    ``ell_terms`` field in the public input is rejected.
    """
    model = data["endpoint_model"]
    if "u_terms" in model or "ell_terms" in model:
        raise ValueError("hand-written u_terms/ell_terms are not accepted; use the endpoint authority")
    convention_path = input_path.parent / model["endpoint_convention"]
    if not convention_path.is_file():
        raise ValueError("missing endpoint convention file")
    convention_sha = sha256_bytes(convention_path.read_bytes())
    if convention_sha != model["endpoint_convention_sha256"]:
        raise ValueError("endpoint convention SHA mismatch")
    builder = HERE / "build_t73_endpoint_transport.py"
    spec = importlib.util.spec_from_file_location("build_t73_endpoint_transport", builder)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load build_t73_endpoint_transport.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    derived = module.derive_endpoint_terms(convention_path)
    derived["convention_file_sha256"] = convention_sha
    return derived


def build_receipt(input_path: Path) -> dict[str, Any]:
    input_bytes = input_path.read_bytes()
    data = json.loads(input_bytes.decode("utf-8"))
    if data.get("schema") != "t73_delta3_public_input/v1":
        raise ValueError("unsupported input schema")
    check_binding(data)

    degree = data["endpoint_model"]["truncation_degree"]
    dimension = data["endpoint_model"]["dimension"]
    if (degree, dimension) != (6, 88):
        raise ValueError("unexpected endpoint model")
    if data["endpoint_model"]["burau_positive_block"] != [["1-t", "t"], ["1", "0"]]:
        raise ValueError("unexpected Burau convention")
    if data["endpoint_model"]["parameter_substitution"] != "t=q^-2; h=q-1":
        raise ValueError("unexpected parameter substitution")
    position_table_path = input_path.parent / "B88_POSITION_TO_PASSAGE_TABLE.json"
    if not position_table_path.is_file():
        raise ValueError("missing B88 position-to-passage table")
    position_table_sha = sha256_bytes(position_table_path.read_bytes())
    if position_table_sha != data["endpoint_model"]["position_table_sha256"]:
        raise ValueError("B88 position-to-passage table SHA mismatch")

    b44, _ = build_oriented_b44(data)
    b88 = cable_word(b44)
    integrity = data["point_push"]["derived_integrity"]
    if len(b88) != integrity["B88_length"]:
        raise ValueError("B88 length mismatch")
    if canonical_sha(b88) != integrity["B88_sha256"]:
        raise ValueError("B88 SHA mismatch")

    endpoint = derive_endpoint_terms_from_authority(input_path, data)
    u_terms = endpoint["u_terms"]
    ell_terms = endpoint["ell_terms"]
    vector = sparse_vector(dimension, degree, u_terms)
    delta_vector = delta_apply(b88, vector)
    eta_epsilon = apply_covector(delta_vector, ell_terms)
    delta_xi_vector = delta_apply(b88, delta_vector)
    xi_epsilon = apply_covector(delta_xi_vector, ell_terms)
    eta_h = substitute_epsilon_with_h(eta_epsilon, degree)
    xi_h = substitute_epsilon_with_h(xi_epsilon, degree)

    return {
        "schema": "t73_delta3_public_receipt/v2",
        "input_sha256": sha256_bytes(input_bytes),
        "endpoint_model": {
            "u_terms": u_terms,
            "ell_terms": ell_terms,
            "derived_from": (
                "data/T73_ENDPOINT_CONVENTION.json via "
                "scripts/build_t73_endpoint_transport.py derive_endpoint_terms"
            ),
            "endpoint_convention_sha256": endpoint["convention_file_sha256"],
            "position_table_sha256": position_table_sha,
        },
        "derived_words": {
            "B44_length": len(b44),
            "B44_sha256": canonical_sha(b44),
            "B88_length": len(b88),
            "B88_sha256": canonical_sha(b88),
        },
        "epsilon_of_h_degrees_0_to_6": epsilon_of_h(degree),
        "exact_epsilon_scalars": {
            "ell_(rhoW-I)_u_degrees_0_to_6": eta_epsilon,
            "ell_(rhoW-I)^2_u_degrees_0_to_6": xi_epsilon,
        },
        "exact_h_scalars": {
            "ell_(rhoW-I)_u_degrees_0_to_6": eta_h,
            "ell_(rhoW-I)^2_u_degrees_0_to_6": xi_h,
        },
        "results": {
            "delta3_eta_R_T1": eta_h[3],
            "delta3_xi": xi_h[3],
            "plain_shadow_cubic_of_xi": eta_h[3],
        },
    }


def print_text(receipt: dict[str, Any]) -> None:
    words = receipt["derived_words"]
    eps = receipt["exact_epsilon_scalars"]
    results = receipt["results"]
    print(f"INPUT_SHA256={receipt['input_sha256']}")
    print(f"B44_LENGTH={words['B44_length']}")
    print(f"B44_SHA256={words['B44_sha256']}")
    print(f"B88_LENGTH={words['B88_length']}")
    print(f"B88_SHA256={words['B88_sha256']}")
    print(f"POSITION_TABLE_SHA256={receipt['endpoint_model']['position_table_sha256']}")
    print(f"ENDPOINT_CONVENTION_SHA256={receipt['endpoint_model']['endpoint_convention_sha256']}")
    print("U_TERMS_DERIVED=" + json.dumps(receipt["endpoint_model"]["u_terms"], separators=(",", ":")))
    print("ELL_TERMS_DERIVED=" + json.dumps(receipt["endpoint_model"]["ell_terms"], separators=(",", ":")))
    print(
        "ELL_RHOW_MINUS_I_U_EPS="
        + json.dumps(eps["ell_(rhoW-I)_u_degrees_0_to_6"], separators=(",", ":"))
    )
    print(
        "ELL_RHOW_MINUS_I_SQUARED_U_EPS="
        + json.dumps(
            eps["ell_(rhoW-I)^2_u_degrees_0_to_6"], separators=(",", ":")
        )
    )
    print(f"DELTA3_ETA_T1={results['delta3_eta_R_T1']}")
    print(f"DELTA3_XI={results['delta3_xi']}")
    print(f"PLAIN_SHADOW_CUBIC_XI={results['plain_shadow_cubic_of_xi']}")
    print("VERIFY=PASS")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--json", action="store_true", help="print the receipt as JSON")
    parser.add_argument("--check", action="store_true", help="accepted for CI readability")
    parser.add_argument(
        "--write-receipt",
        type=Path,
        nargs="?",
        const=DEFAULT_RECEIPT,
        help="write the exact JSON receipt (default: data/T73_DELTA3_PUBLIC_RECEIPT.json)",
    )
    args = parser.parse_args()
    receipt = build_receipt(args.input.resolve())
    if args.write_receipt is not None:
        output_path = args.write_receipt.resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
            newline="\n",
        )
    if args.json:
        print(json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True))
    else:
        print_text(receipt)


if __name__ == "__main__":
    main()
