"""
Term-wise decompression of quad-compressed symbol data.

Unlike expand_symb() (which does algebraic basis expansion via front/back space vectors),
UnQuad works prefix-by-prefix, using explicit physics relations (final-entry,
triple-adjacency, integrability) to reconstruct all 4-letter endings from the 8
stored quad basis entries per prefix.
"""

# ── Constants ────────────────────────────────────────────────────────────────

quad_bases = ["dddd", "bbbd", "bdbd", "bbdd", "dbdd", "fbdd", "dbbd", "cddd"]

# Map each basis ending to the repo's @QUAD_i token format (0-indexed)
quad_codes = {quad_bases[i]: f"@QUAD_{i}" for i in range(len(quad_bases))}

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


def UnQuad(prefix, data=None):
    """Decompress all entries with a given (2L-4)-length prefix.

    Args:
        prefix: The symbol prefix of length 2L-4.
        data: Pre-loaded quad data dict (from Phi2Symb(loop, "quad")).
              If None, loads automatically.

    Returns:
        Dict mapping full symbol keys to integer coefficients.
    """
    from aiamplitudes_common_public import Phi2Symb

    res = {}
    loop = len(prefix) // 2 + 2
    if data is None:
        data = Phi2Symb(loop, "quad")

    # Read the 24 base keys associated to this prefix
    linquads, orig = get24(prefix + "dddd")
    for i, d in enumerate(linquads):
        k = d[:-4] + quad_codes[d[-4:]]
        if k in data:
            v = data[k]
            res[orig[i]] = v

    # Apply eq4: 4-letter suffix relations
    for e in eq4:
        key = prefix + e[0]
        val = 0
        for j in range(len(e[1]) // 2):
            k = prefix + e[1][2 * j + 1]
            if k in res:
                val += e[1][2 * j] * res[k]
        if round(val) != 0:
            res[key] = round(val)

    # Apply eq3: 3-letter suffix relations (with each letter prepended)
    for e in eq3:
        for l in ['a', 'b', 'c', 'd', 'e', 'f']:
            key = prefix + l + e[0]
            if Admissible(key):
                val = 0
                for j in range(len(e[1]) // 2):
                    k = prefix + l + e[1][2 * j + 1]
                    if k in res:
                        val += e[1][2 * j] * res[k]
                if round(val) != 0:
                    res[key] = round(val)

    # Apply eq2: 2-letter suffix relations (with each valid 2-letter pair prepended)
    # Bug fix: original notebook had 'ab''ac' (string concat) instead of 'ab','ac'
    eq2_prefixes = [
        'aa', 'ab', 'ac', 'ae', 'af',
        'ba', 'bb', 'bc', 'bd', 'bf',
        'ca', 'cb', 'cc', 'cd', 'ce',
        'db', 'dc', 'dd',
        'ea', 'ec', 'ee',
        'fa', 'fb', 'ff',
    ]
    for e in eq2:
        for l in eq2_prefixes:
            key = prefix + l + e[0]
            if Admissible(key):
                val = 0
                for j in range(len(e[1]) // 2):
                    k = prefix + l + e[1][2 * j + 1]
                    if k in res:
                        val += e[1][2 * j] * res[k]
                if round(val) != 0:
                    res[key] = round(val)
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
