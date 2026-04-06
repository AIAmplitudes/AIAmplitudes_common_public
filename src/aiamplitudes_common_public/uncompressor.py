"""
Term-wise decompression of compressed symbol data.

UnQuad maps 24 quad-basis indep values to all 4-letter suffix values using a
matrix derived (lazily, on first use) from the eq4/eq3/eq2 physics relations.
The matrix depends on the last letter of the prefix (6 variants, cached).

UnOct maps 279 oct-basis indep values (93 bases × 3 rotations) to all
8-letter suffix values using precomputed data (data/oct_matrices.npz).
Two-stage exact integer arithmetic:
  1. C^{-T} converts oct-lookup values to the restrictive B-space basis
  2. Sparse B matrix maps basis values to suffix values (per-column exact division)

Single-term lookups (UnQuadTerm, UnOctTerm) compute a single suffix value
via dot products instead of decompressing all suffixes.
"""

# ── Constants ────────────────────────────────────────────────────────────────

quad_bases = ["dddd", "bbbd", "bdbd", "bbdd", "dbdd", "fbdd", "dbbd", "cddd"]

# Map each basis ending to the repo's @QUAD_i token format (0-indexed)
quad_codes = {b: f"@QUAD_{i}" for i, b in enumerate(quad_bases)}

# 16 dihedral equivalents of the quad bases
quad_cyclic = [
    "eeee", "ffff", "ccce", "aaaf", "cece", "afaf", "ccee", "aaff",
    "ecee", "faff", "dcee", "eaff", "ecce", "faaf", "aeee", "bfff",
]

AdjacencyMap = [
    [1, 1, 1, 0, 1, 1],
    [1, 1, 1, 1, 0, 1],
    [1, 1, 1, 1, 1, 0],
    [0, 1, 1, 1, 0, 0],
    [1, 0, 1, 0, 1, 0],
    [1, 1, 0, 0, 0, 1],
]

# ── Relation tables ──────────────────────────────────────────────────────────
# 18 four-letter suffix relations (final-entry relations at weight 4)
eq4 = [
    ["bddd", [-1, "faff", 1, "dbdd", 1, "eaff", -1, "fbdd", 1, "aeee"]],
    ["ceee", [-1, "dbdd", 1, "ecee", 1, "fbdd", -1, "dcee", 1, "bfff"]],
    ["afff", [-1, "ecee", 1, "faff", 1, "dcee", -1, "eaff", 1, "cddd"]],
    ["abdd", [-0.5, "ecee", 0.5, "faff", -0.5, "eaff", 0.5, "dcee", 0.5, "cddd", -0.5, "aeee"]],
    ["bcee", [-0.5, "fbdd", 0.5, "dbdd", -0.5, "bfff", -0.5, "faff", 0.5, "eaff", 0.5, "aeee"]],
    ["caff", [-0.5, "dcee", 0.5, "ecee", -0.5, "cddd", -0.5, "dbdd", 0.5, "fbdd", 0.5, "bfff"]],
    ["cbdd", [-0.5, "dcee", 0.5, "ecee", -0.5, "cddd", -0.5, "dbdd", 0.5, "fbdd", 0.5, "bfff"]],
    ["acee", [-0.5, "ecee", 0.5, "faff", -0.5, "eaff", 0.5, "dcee", 0.5, "cddd", -0.5, "aeee"]],
    ["baff", [-0.5, "fbdd", 0.5, "dbdd", -0.5, "bfff", -0.5, "faff", 0.5, "eaff", 0.5, "aeee"]],
    ["cdbd", [-0.5, "dcee", 0.5, "ecee", -0.5, "cddd", -0.5, "dbdd", 0.5, "fbdd", 0.5, "bfff"]],
    ["aece", [-0.5, "ecee", 0.5, "faff", -0.5, "eaff", 0.5, "dcee", 0.5, "cddd", -0.5, "aeee"]],
    ["bfaf", [-0.5, "fbdd", 0.5, "dbdd", -0.5, "bfff", -0.5, "faff", 0.5, "eaff", 0.5, "aeee"]],
    ["ddbd", [1, "dbdd"]],
    ["eece", [1, "ecee"]],
    ["ffaf", [1, "faff"]],
    ["fbbd", [-1, "bbdd", 1, "dbbd", 0.5, "faff", -0.5, "dbdd", 0.5, "fbdd", -0.5, "eaff", -0.5, "aeee", 0.5, "bfff"]],
    ["dcce", [-1, "ccee", 1, "ecce", 0.5, "dbdd", -0.5, "ecee", 0.5, "dcee", -0.5, "fbdd", -0.5, "bfff", 0.5, "cddd"]],
    ["eaaf", [-1, "aaff", 1, "faaf", 0.5, "ecee", -0.5, "faff", 0.5, "eaff", -0.5, "dcee", -0.5, "cddd", 0.5, "aeee"]],
]

# 6 three-letter suffix relations (applied with each letter prepended)
eq3 = [
    ["cdd", [-1, "cee"]],
    ["bff", [-1, "bdd"]],
    ["aee", [-1, "aff"]],
    ["fbd", [1, "dbd", -1, "bdd"]],
    ["dce", [1, "ece", -1, "cee"]],
    ["eaf", [1, "faf", -1, "aff"]],
]

# 3 two-letter suffix relations (applied with each valid 2-letter pair prepended)
eq2 = [
    ["bf", [1, "bd"]],
    ["ae", [1, "af"]],
    ["cd", [1, "ce"]],
]


# ── Utility functions ────────────────────────────────────────────────────────

def Subst(key, transf):
    """Apply a 6-letter permutation (string transf) to key."""
    res = ""
    for l in key:
        res += transf[ord(l) - ord('a')]
    return res


def DihedralEq(key):
    """Return the set of (up to 6) dihedral equivalents of key."""
    result = set()
    rotations = ["cabfde", "bcaefd", "abcdef"]
    flip = "bacedf"
    for i in range(3):
        result.add(Subst(key, rotations[i]))
        result.add(Subst(Subst(key, rotations[i]), flip))
    return result


def Admissible(key):
    """Check that key satisfies initial, final, and adjacency conditions."""
    if len(key) < 2 or len(key) % 2 != 0:
        return False
    if key[0] in ['d', 'e', 'f'] or key[-1] in ['a', 'b', 'c']:
        return False
    for i in range(len(key) - 1):
        if AdjacencyMap[ord(key[i]) - ord('a')][ord(key[i + 1]) - ord('a')] == 0:
            return False
    return True


def UnDihedral(key, ends="bf"):
    """Return the canonical dihedral representative (ending in af*)."""
    rotations = ["cabfde", "bcaefd", "abcdef"]
    flip = "bacedf"
    rotkey = Subst(key, rotations[ord(key[-1]) - ord('d')])
    pos = len(rotkey) - 1
    while rotkey[pos] == 'f':
        pos -= 1
    if rotkey[pos] == 'b':
        return Subst(rotkey, flip)
    return rotkey


# ── Core functions ───────────────────────────────────────────────────────────

def get24(key):
    """Given a key (prefix + 4-letter ending), return the 24 quad-basis keys
    and their dihedral-equivalent originals.

    Returns (linquads, orig) where linquads[i] has the basis ending suitable
    for lookup in the quad file, and orig[i] is the actual symbol key.
    """
    res = []
    res2 = []
    rotations = ["cabfde", "bcaefd"]
    for i in range(len(quad_bases)):
        res.append(key[:-4] + quad_bases[i])
        res2.append(key[:-4] + quad_bases[i])
        res.append(Subst(key[:-4] + quad_cyclic[2 * i], rotations[0]))
        res2.append(key[:-4] + quad_cyclic[2 * i])
        res.append(Subst(key[:-4] + quad_cyclic[2 * i + 1], rotations[1]))
        res2.append(key[:-4] + quad_cyclic[2 * i + 1])
    return res, res2


# Valid 2-letter prefixes for eq2 relations
eq2_prefixes = [
    'aa', 'ab', 'ac', 'ae', 'af',
    'ba', 'bb', 'bc', 'bd', 'bf',
    'ca', 'cb', 'cc', 'cd', 'ce',
    'db', 'dc', 'dd',
    'ea', 'ec', 'ee',
    'fa', 'fb', 'ff',
]

# 24 indep suffixes in get24 order: for each quad basis i, the canonical
# suffix followed by its two dihedral copies
_indep_order = [s for j in range(8)
                for s in (quad_bases[j], quad_cyclic[2*j], quad_cyclic[2*j+1])]

# Minimal valid prefixes for each possible last letter (for matrix derivation)
_min_prefix = {'a': 'aa', 'b': 'ab', 'c': 'ac', 'd': 'bd', 'e': 'ae', 'f': 'af'}

_quad_matrix_cache = {}

def _apply_relations_exact(prefix, res):
    """Apply eq4, eq3, eq2 relations.

    Used for matrix derivation where indep values are 0/1 and intermediate
    values can be half-integers (due to ±0.5 coefficients in eq4).
    """
    for e in eq4:
        key = prefix + e[0]
        val = 0
        for j in range(len(e[1]) // 2):
            k = prefix + e[1][2 * j + 1]
            if k in res:
                val += e[1][2 * j] * res[k]
        if val != 0:
            res[key] = val

    for e in eq3:
        for l in 'abcdef':
            key = prefix + l + e[0]
            if Admissible(key):
                val = 0
                for j in range(len(e[1]) // 2):
                    k = prefix + l + e[1][2 * j + 1]
                    if k in res:
                        val += e[1][2 * j] * res[k]
                if val != 0:
                    res[key] = val

    for e in eq2:
        for l in eq2_prefixes:
            key = prefix + l + e[0]
            if Admissible(key):
                val = 0
                for j in range(len(e[1]) // 2):
                    k = prefix + l + e[1][2 * j + 1]
                    if k in res:
                        val += e[1][2 * j] * res[k]
                if val != 0:
                    res[key] = val


def _get_quad_matrix(last_letter):
    """Return (B, suffixes, sf_index) for the given last prefix letter.

    B is (24, N_suffixes) such that indep_values @ B gives all suffix values.
    Cached after first computation.
    """
    if last_letter in _quad_matrix_cache:
        return _quad_matrix_cache[last_letter]

    import numpy as np

    prefix = _min_prefix[last_letter]
    plen = len(prefix)

    # Probe with each unit vector to build the matrix
    all_suffixes = set()
    columns = []
    for sf in _indep_order:
        res = {prefix + sf: 1}
        _apply_relations_exact(prefix, res)
        col = {k[plen:]: v for k, v in res.items()}
        columns.append(col)
        all_suffixes.update(col.keys())

    suffixes = sorted(all_suffixes)
    sf_index = {sf: j for j, sf in enumerate(suffixes)}

    B = np.zeros((24, len(suffixes)))
    for i, col in enumerate(columns):
        for sf, v in col.items():
            B[i, sf_index[sf]] = v

    _quad_matrix_cache[last_letter] = (B, suffixes, sf_index)
    return B, suffixes, sf_index


def UnQuad(prefix, data=None):
    """Decompress all entries with a given (2L-4)-length prefix.

    Uses a precomputed matrix (derived from the eq4/eq3/eq2 relations)
    to map the 24 quad-basis indep values to all suffix values in one
    matrix multiply.

    Args:
        prefix: The symbol prefix of length 2L-4.
        data: Pre-loaded quad data dict (from Phi2Symb(loop, "quad")).
              If None, loads automatically.

    Returns:
        Dict mapping full symbol keys to integer coefficients.
    """
    import numpy as np
    from aiamplitudes_common_public import Phi2Symb

    loop = len(prefix) // 2 + 2
    if data is None:
        data = Phi2Symb(loop, "quad")

    # Read the 24 indep values via get24
    linquads, orig = get24(prefix + "dddd")
    indep_values = np.zeros(24)
    for i, d in enumerate(linquads):
        k = d[:-4] + quad_codes[d[-4:]]
        if k in data:
            indep_values[i] = data[k]

    # Matrix multiply: indeps -> all suffixes
    B, suffixes, _ = _get_quad_matrix(prefix[-1])
    vals = indep_values @ B

    res = {}
    for j, sf in enumerate(suffixes):
        v = round(vals[j])
        if v != 0:
            res[prefix + sf] = v
    return res


def UnQuadLoop(loop):
    """Decompress an entire loop order from quad format.

    Loads the quad data once, collects all prefixes and their dihedral
    equivalents, then calls UnQuad on each.

    Args:
        loop: Loop number (integer).

    Returns:
        Dict mapping all full symbol keys to integer coefficients.
    """
    from aiamplitudes_common_public import Phi2Symb

    res = {}
    data = Phi2Symb(loop, "quad")
    todo = set()
    # Collect all prefixes and their dihedral equivalents
    for d in data:
        # Strip the @QUAD_i token (find @ and take everything before it)
        prefix = d[:d.index('@')]
        todo.update(DihedralEq(prefix))
    # UnQuad each prefix
    for t in todo:
        if len(t) != 2 * loop - 4:
            print("aie", t)
        res.update(UnQuad(t, data=data))
    return res

# ── Single-term lookup ──────────────────────────────────────────────────────

def UnQuadTerm(key, data=None):
    """Look up a single symbol term from quad-compressed data.

    Instead of decompressing all suffixes for the prefix, computes just
    the dot product of the 24 indep values with the relevant column of the
    quad basis matrix.

    Args:
        key: Full symbol key (length 2L).
        data: Pre-loaded quad data dict. If None, loads automatically.

    Returns:
        Integer coefficient of this term (0 if absent).
    """
    import numpy as np
    from aiamplitudes_common_public import Phi2Symb

    loop = len(key) // 2
    prefix = key[:-4]
    suffix = key[-4:]

    if data is None:
        data = Phi2Symb(loop, "quad")

    B, suffixes, sf_index = _get_quad_matrix(prefix[-1])
    if suffix not in sf_index:
        return 0
    j = sf_index[suffix]

    # Read the 24 indep values
    linquads, _ = get24(prefix + "dddd")
    indep_values = np.zeros(24)
    for i, d in enumerate(linquads):
        k = d[:-4] + quad_codes[d[-4:]]
        if k in data:
            indep_values[i] = data[k]

    return round(indep_values @ B[:, j])


# ── Oct (octuple) decompression ────────────────────────────────────────────

oct_bases = [
    "aaaaaaaf", "aaaaafaf", "aaaafaaf", "aaafaaaf", "aaafaaff", "aaafafaf",
    "aaeaffff", "aaeeaaaf", "aaeeaaff", "aaeeafaf", "aafaaaaf", "aafaafaf",
    "aafafaaf", "aafaffff", "aafbbdbd", "aafbdddd", "aafbffff", "aaffbdbd",
    "aaffffff", "aeaaaaaf", "aeaaafaf", "aeafaaaf", "aeafafaf", "aeeeaaaf",
    "aeeeaaff", "aeeecece", "afaaaaaf", "afaaafaf", "afaafaaf", "affaaaaf",
    "affaafaf", "affafaaf", "affbbdbd", "affbdbbd", "afffaaaf", "afffbdbd",
    "bafbffff", "bbaaffff", "bbaeaaaf", "bbbabdbd", "bbbaeeee", "bbbbcece",
    "bbbddbdd", "bbbfbbbd", "bbbfbbdd", "bbbfbdbd", "bbcbdddd", "bbccecce",
    "bbddbbbd", "bbfaafaf", "bbfbbdbd", "bbfbdbbd", "bbfbffff", "bbffffff",
    "bcdbdddd", "bcddbbdd", "bceeeeee", "bdbabdbd", "bdbaeeee", "bdbddbdd",
    "bddbdddd", "bddddddd", "dbbbbbbd", "dbbbbdbd", "dbbbdbbd", "dbbdbbdd",
    "dbbddbdd", "dbbfbbbd", "dbbfbdbd", "dbccecce", "dbcceeee", "dbdbdddd",
    "dbdcdddd", "dbdddddd", "dbfaffff", "dbfbbbbd", "dbfbbdbd", "dbfbdbbd",
    "dcaeafaf", "dcaeeeee", "dcafafaf", "dcbaeeee", "dccceeee", "dccdccee",
    "dccecece", "dcdcdddd", "dcddcece", "dcecdddd", "dceeccce", "dceeccee",
    "dceeeeee", "ddbbbdbd", "dddddddd",
]

oct_codes = {b: f"@OCT_{i}" for i, b in enumerate(oct_bases)}

# Cyclic images of oct_bases: oct_cyclic[2*i] = rot2(base_i),
# oct_cyclic[2*i+1] = rot1(base_i).  This ordering matches quad_cyclic
# so that get279 can follow the same pattern as get24: applying rot1 to
# (prefix + cyclic[2*i]) recovers base_i in the last 8 chars.
oct_cyclic = [s for ob in oct_bases
              for s in (Subst(ob, "bcaefd"), Subst(ob, "cabfde"))]  # rot2, rot1

# ── Oct lookup and decompression ────────────────────────────────────────


def get279(prefix):
    """Return the 279 oct-basis lookup keys for a given prefix.

    For each of 93 oct bases, returns the key to look up in the oct data
    file for the identity rotation and the two non-identity rotations
    (which access the rotated prefix's data).

    Args:
        prefix: The symbol prefix of length 2L-8.

    Returns:
        (linkeys, orig) where linkeys[i] is the key to look up in the oct
        data (with @OCT_j suffix), and orig[i] is the actual full symbol key.
    """
    rotations = ["cabfde", "bcaefd"]
    linkeys = []
    orig = []
    for i in range(len(oct_bases)):
        # Identity: look up this prefix's oct data for base i
        linkeys.append(prefix + oct_codes[oct_bases[i]])
        orig.append(prefix + oct_bases[i])
        # Rotation 1: rotate(prefix + cyclic_2i) -> look up rotated prefix's data
        rot1_word = Subst(prefix + oct_cyclic[2 * i], rotations[0])
        linkeys.append(rot1_word[:-8] + oct_codes[rot1_word[-8:]])
        orig.append(prefix + oct_cyclic[2 * i])
        # Rotation 2: rotate(prefix + cyclic_2i+1) -> look up rotated prefix's data
        rot2_word = Subst(prefix + oct_cyclic[2 * i + 1], rotations[1])
        linkeys.append(rot2_word[:-8] + oct_codes[rot2_word[-8:]])
        orig.append(prefix + oct_cyclic[2 * i + 1])
    return linkeys, orig


_oct_matrix_cache = {}


def _load_oct_matrices():
    """Load precomputed oct decompression data from data/oct_matrices.npz.

    Two-stage decompression:
      1. rest_vals = C_inv_T @ oct_lookup  (279→279, exact rational→integer)
      2. suffix_values = rest_vals @ B_sparse  (per-column exact division)
    """
    if _oct_matrix_cache:
        return
    import numpy as np
    import os
    datadir = os.path.join(os.path.dirname(__file__), '..', '..', 'data')
    path = os.path.join(datadir, 'oct_matrices.npz')
    mdata = np.load(path)

    # C_inv_T: (279×279) stored as per-row scaled integers
    _oct_matrix_cache['C_inv_T_num'] = mdata['C_inv_T_num']      # (279,279) int64
    _oct_matrix_cache['C_inv_T_denoms'] = mdata['C_inv_T_denoms']  # (279,) int64

    for ll in 'abcdef':
        suffixes = list(mdata[f'suffixes_{ll}'])
        sf_index = {sf: j for j, sf in enumerate(suffixes)}
        _oct_matrix_cache[ll] = {
            'suffixes': suffixes,
            'sf_index': sf_index,
            'rows': mdata[f'rows_{ll}'],        # int32, CSC row indices
            'nums': mdata[f'nums_{ll}'],        # int64, scaled numerators
            'col_ptrs': mdata[f'col_ptrs_{ll}'],  # int32, column pointers
            'col_lcms': mdata[f'col_lcms_{ll}'],  # int64, per-column denominators
        }


def _oct_step1(oct_lookup):
    """Step 1: Convert oct_lookup (279 ints) to rest_vals (279 ints) via C_inv_T.

    Uses exact integer arithmetic with per-row denominators.
    """
    import numpy as np
    _load_oct_matrices()
    C_num = _oct_matrix_cache['C_inv_T_num']
    C_den = _oct_matrix_cache['C_inv_T_denoms']
    oct_arr = np.asarray(oct_lookup, dtype=np.int64)
    raw = C_num @ oct_arr
    assert np.all(raw % C_den == 0), "Non-integer rest_vals in C_inv_T step"
    return (raw // C_den).tolist()


def _oct_step2(prefix, rest_vals):
    """Step 2: Multiply rest_vals (279 ints) by sparse B to get suffix values.

    Uses per-column common denominators for exact division.
    """
    import numpy as np
    _load_oct_matrices()
    m = _oct_matrix_cache[prefix[-1]]
    suffixes = m['suffixes']
    rows = m['rows']
    nums = m['nums']
    col_ptrs = m['col_ptrs']
    col_lcms = m['col_lcms']
    n_suf = len(suffixes)
    rest_arr = np.asarray(rest_vals, dtype=np.int64)

    res = {}
    for j in range(n_suf):
        start, end = int(col_ptrs[j]), int(col_ptrs[j + 1])
        if start == end:
            continue
        total = int(np.dot(rest_arr[rows[start:end]], nums[start:end]))
        col_lcm = int(col_lcms[j])
        if col_lcm != 1:
            assert total % col_lcm == 0, \
                f"Non-integer at suffix {suffixes[j]}: {total}/{col_lcm}"
            total //= col_lcm
        if total != 0:
            res[prefix + suffixes[j]] = total
    return res


def UnOct(prefix, data=None):
    """Decompress all entries with a given (2L-8)-length prefix.

    Two-stage exact decompression:
      1. Convert oct-lookup values to restrictive-basis values via C^{-T}
      2. Multiply by sparse B matrix with per-column exact division

    Args:
        prefix: The symbol prefix of length 2L-8.
        data: Pre-loaded oct data dict (from Phi2Symb(loop, "oct")).
    Returns:
        Dict mapping full symbol keys to integer coefficients.
    """
    from aiamplitudes_common_public import Phi2Symb

    loop = len(prefix) // 2 + 4
    if data is None:
        data = Phi2Symb(loop, "oct")

    # Read the 279 oct-lookup values via get279
    linkeys, _ = get279(prefix)
    oct_lookup = [0] * 279
    for i, k in enumerate(linkeys):
        if k in data:
            oct_lookup[i] = data[k]

    # Step 1: oct_lookup -> rest_vals
    rest_vals = _oct_step1(oct_lookup)

    # Step 2: rest_vals -> suffix values
    return _oct_step2(prefix, rest_vals)


def UnOctLoop(loop):
    """Decompress an entire loop order from oct format.

    Loads the oct data once, collects all prefixes and their dihedral
    equivalents, then calls UnOct on each.

    Args:
        loop: Loop number (integer).

    Returns:
        Dict mapping all full symbol keys to integer coefficients.
    """
    from aiamplitudes_common_public import Phi2Symb

    res = {}
    data = Phi2Symb(loop, "oct")
    todo = set()
    for d in data:
        prefix = d[:d.index('@')]
        todo.update(DihedralEq(prefix))
    for t in todo:
        if len(t) != 2 * loop - 8:
            print("aie", t)
        res.update(UnOct(t, data=data))
    return res


def UnOctTerm(key, data=None):
    """Look up a single symbol term from oct-compressed data.

    Uses the same two-stage exact decompression as UnOct but only
    computes the single requested suffix.

    Args:
        key: Full symbol key (length 2L).
        data: Pre-loaded oct data dict. If None, loads automatically.

    Returns:
        Integer coefficient of this term (0 if absent).
    """
    from aiamplitudes_common_public import Phi2Symb

    loop = len(key) // 2
    prefix = key[:-8]
    suffix = key[-8:]

    if data is None:
        data = Phi2Symb(loop, "oct")

    _load_oct_matrices()
    m = _oct_matrix_cache[prefix[-1]]
    if suffix not in m['sf_index']:
        return 0
    j = m['sf_index'][suffix]

    # Read oct-lookup values
    linkeys, _ = get279(prefix)
    oct_lookup = [0] * 279
    for i, k in enumerate(linkeys):
        if k in data:
            oct_lookup[i] = data[k]

    # Step 1: oct_lookup -> rest_vals
    rest_vals = _oct_step1(oct_lookup)

    # Step 2: single column dot product
    rows = m['rows']
    nums = m['nums']
    col_ptrs = m['col_ptrs']
    start = int(col_ptrs[j])
    end = int(col_ptrs[j + 1])
    col_lcm = int(m['col_lcms'][j])

    total = 0
    for idx in range(start, end):
        total += rest_vals[int(rows[idx])] * int(nums[idx])
    if col_lcm != 1:
        assert total % col_lcm == 0
        total //= col_lcm
    return total
