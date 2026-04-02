"""
Derive the OCT decompression data analytically from B-space relations.

Produces two components stored in data/oct_matrices.npz:

1. C_inv_T (279×279 exact rational): converts oct_lookup values (93 × 3 rotations)
   to the restrictive B-space basis values. Stored as int numerators + per-row denoms.

2. B_rest (sparse, per prefix-last-letter): maps restrictive-basis values to all
   8-letter suffix values. Stored as sparse (row, col, numerator) with per-column
   common denominators.

At decompression time:
  rest_vals = C_inv_T @ oct_lookup          (exact integer result)
  suffix_values[j] = sum_i rest_vals[i] * B[i,j]  (per-column exact division)
"""
import sys
import os
import time
import re
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from aiamplitudes_common_public.uncompressor import (
    _load_oct_rels, _get_prepends, _min_prefix, Subst, Admissible
)
from aiamplitudes_common_public.fbspaces import get_rest_bspace, get_brels
from aiamplitudes_common_public.download_data import relpath
from fractions import Fraction
from math import gcd
from collections import defaultdict


def lcm(a, b):
    return a * b // gcd(a, b)


def apply_backspace_only(prefix, res):
    """Back-space relation cascade (weights 8 down to 2)."""
    nonzero = _load_oct_rels()
    plen = len(prefix)
    last_idx = ord(prefix[-1]) - ord('a')
    for dep, combo in nonzero[8].items():
        key = prefix + dep
        if key not in res:
            val = Fraction(0)
            for ind, coef in combo.items():
                k = prefix + ind
                if k in res:
                    val += coef * res[k]
            if val != 0:
                res[key] = val
    for w in range(7, 1, -1):
        prepends = _get_prepends(8 - w, last_idx)
        for dep, combo in nonzero[w].items():
            for pre in prepends:
                key = prefix + pre + dep
                if len(key) != plen + 8 or key in res:
                    continue
                if not Admissible(key):
                    continue
                val = Fraction(0)
                for ind, coef in combo.items():
                    k = prefix + pre + ind
                    if k in res:
                        val += coef * res[k]
                if val != 0:
                    res[key] = val


# ── Parse octic93 ────────────────────────────────────────────────────────

datapath = os.path.join(os.path.dirname(__file__), '..', '..', 'LanceData_New',
                        'oldFiles', 'EZ_symb_oct_new_norm')
with open(datapath, 'r') as f:
    content = f.read()

match = re.search(r'octic93 := \[(.*?)\] :', content, re.DOTALL)
entries = re.findall(r'E\(([^)]+)\)', match.group(1))
octic93 = [''.join(e.replace(',', '').replace(' ', '')) for e in entries]
print(f"Parsed {len(octic93)} oct bases")

# Order matches get279: [base, rot2(base), rot1(base)] per group
oct_expanded = []
for base in octic93:
    oct_expanded.append(base)
    oct_expanded.append(Subst(base, 'bcaefd'))   # rot2
    oct_expanded.append(Subst(base, 'cabfde'))   # rot1

# ── Load basis and relations ─────────────────────────────────────────────

flip8, _ = get_rest_bspace(8)
rest_basis = [v for k, v in sorted(flip8.items(), key=lambda x: int(x[0].split('_')[2]))]
rest_index = {w: i for i, w in enumerate(rest_basis)}
rest_set = set(rest_basis)
n = 279

all_nonzero = {}
all_zero = {}
for w in range(2, 9):
    rels = get_brels(w, relpath)
    all_nonzero[w] = {k: v for k, v in rels.items() if v != {None: 0}}
    all_zero[w] = {k for k, v in rels.items() if v == {None: 0}}


def express_in_basis(word):
    if word in rest_set:
        return {word: Fraction(1)}
    if word in all_nonzero[8]:
        return dict(all_nonzero[8][word])
    if word in all_zero[8]:
        return {}
    for w in range(7, 1, -1):
        dep = word[8 - w:]
        pre = word[:8 - w]
        if dep in all_nonzero[w]:
            combo = all_nonzero[w][dep]
            result = {}
            for ind, coef in combo.items():
                sub = express_in_basis(pre + ind)
                for bw, bc in sub.items():
                    result[bw] = result.get(bw, Fraction(0)) + coef * bc
            return result
        if dep in all_zero[w]:
            return {}
    return None


# ── Build and invert C^T ─────────────────────────────────────────────────

print("\nBuilding C^T...")
t0 = time.time()
CT = [[Fraction(0)] * n for _ in range(n)]
for j, word in enumerate(oct_expanded):
    expansion = express_in_basis(word)
    if expansion:
        for bw, coef in expansion.items():
            if bw in rest_index:
                CT[j][rest_index[bw]] = coef
print(f"  Built in {time.time() - t0:.1f}s")

print("Inverting C^T (Gauss-Jordan, exact)...")
t0 = time.time()
aug = [CT[j][:] + [Fraction(1) if i == j else Fraction(0) for i in range(n)]
       for j in range(n)]
for col in range(n):
    pivot = None
    for row in range(col, n):
        if aug[row][col] != 0:
            pivot = row
            break
    assert pivot is not None, f"Singular at col {col}"
    if pivot != col:
        aug[col], aug[pivot] = aug[pivot], aug[col]
    pv = aug[col][col]
    if pv != 1:
        for k in range(2 * n):
            aug[col][k] /= pv
    for row in range(n):
        if row == col or aug[row][col] == 0:
            continue
        factor = aug[row][col]
        for k in range(2 * n):
            aug[row][k] -= factor * aug[col][k]
    if (col + 1) % 50 == 0:
        print(f"  col {col + 1}/{n} ({time.time() - t0:.1f}s)")

C_inv_T = [[aug[j][n + i] for i in range(n)] for j in range(n)]
print(f"  Inverted in {time.time() - t0:.1f}s")

# Store as per-row scaled integers
row_lcms = []
for i in range(n):
    d = 1
    for j in range(n):
        d = lcm(d, C_inv_T[i][j].denominator)
    row_lcms.append(d)

C_inv_T_num = np.zeros((n, n), dtype=np.int64)
C_inv_T_denoms = np.array(row_lcms, dtype=np.int64)
for i in range(n):
    for j in range(n):
        C_inv_T_num[i, j] = int(C_inv_T[i][j] * row_lcms[i])

print(f"  C_inv_T: max denom={max(row_lcms)}, "
      f"frac rows={sum(1 for d in row_lcms if d > 1)}")

# ── Build sparse B matrices ─────────────────────────────────────────────

results = {}
for last_letter in 'abcdef':
    prefix = _min_prefix[last_letter]
    plen = len(prefix)
    print(f"\n=== B_{last_letter} (prefix='{prefix}') ===")

    t0 = time.time()
    all_suffixes = set()
    sparse_data = []  # (row_i, suffix_str, Fraction_value)
    for i, basis_sf in enumerate(rest_basis):
        full = prefix + basis_sf
        if not Admissible(full):
            continue
        res = {full: Fraction(1)}
        apply_backspace_only(prefix, res)
        for k, v in res.items():
            sf = k[plen:]
            sparse_data.append((i, sf, v))
            all_suffixes.add(sf)
        if (i + 1) % 50 == 0:
            print(f"  {i + 1}/279 ({time.time() - t0:.1f}s)")

    suffixes = sorted(all_suffixes)
    sf_index = {sf: j for j, sf in enumerate(suffixes)}
    n_suf = len(suffixes)

    # Group by column and compute per-column LCMs
    col_entries = defaultdict(list)
    for row, sf, val in sparse_data:
        col_entries[sf_index[sf]].append((row, val))

    col_lcms = np.ones(n_suf, dtype=np.int64)
    for j, entries in col_entries.items():
        d = 1
        for _, v in entries:
            d = lcm(d, v.denominator)
        col_lcms[j] = d

    # Store as CSC-like sparse: arrays of (row, scaled_numerator) per column
    # Use flat arrays with column pointers
    rows_flat = []
    nums_flat = []
    col_ptrs = [0]
    for j in range(n_suf):
        entries = col_entries.get(j, [])
        for row, val in entries:
            rows_flat.append(row)
            nums_flat.append(int(val * int(col_lcms[j])))
        col_ptrs.append(len(rows_flat))

    rows_arr = np.array(rows_flat, dtype=np.int32)
    nums_arr = np.array(nums_flat, dtype=np.int64)
    ptrs_arr = np.array(col_ptrs, dtype=np.int32)

    dt = time.time() - t0
    int_cols = np.sum(col_lcms == 1)
    print(f"  {n_suf} suffixes, {len(sparse_data)} nnz ({len(sparse_data)/(n*n_suf)*100:.1f}%), "
          f"int cols={int_cols}/{n_suf}, {dt:.1f}s")

    results[last_letter] = {
        'suffixes': suffixes,
        'rows': rows_arr,
        'nums': nums_arr,
        'col_ptrs': ptrs_arr,
        'col_lcms': col_lcms,
    }

# ── Save ─────────────────────────────────────────────────────────────────

outpath = os.path.join(os.path.dirname(__file__), '..', 'data', 'oct_matrices.npz')
os.makedirs(os.path.dirname(outpath), exist_ok=True)

save_dict = {
    'octic93': np.array(octic93),
    'C_inv_T_num': C_inv_T_num,
    'C_inv_T_denoms': C_inv_T_denoms,
}
for ll in 'abcdef':
    r = results[ll]
    save_dict[f'suffixes_{ll}'] = np.array(r['suffixes'])
    save_dict[f'rows_{ll}'] = r['rows']
    save_dict[f'nums_{ll}'] = r['nums']
    save_dict[f'col_ptrs_{ll}'] = r['col_ptrs']
    save_dict[f'col_lcms_{ll}'] = r['col_lcms']

np.savez_compressed(outpath, **save_dict)
filesize = os.path.getsize(outpath)
print(f"\nSaved to {outpath} ({filesize / 1024 / 1024:.1f} MB)")
