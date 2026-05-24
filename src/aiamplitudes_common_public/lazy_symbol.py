"""
Lazy symbol interface for high-loop data.

Provides dict-like access to symbol coefficients without holding the full
uncompressed symbol in memory. Backed by termwise decompression with a
per-prefix cache of the expensive intermediate (rest_vals for oct,
indep_values for quad) so repeated lookups against keys sharing a prefix
skip the matrix work and only pay the cheap per-suffix column dot product.
"""

from aiamplitudes_common_public.uncompressor import (
    _quad_indep_values, _quad_term_from_indep,
    _oct_rest_vals, _oct_term_from_rest_vals,
    _load_oct_matrices, _oct_matrix_cache,
    _get_quad_matrix,
)


class LazySymbol:
    """Dict-like wrapper around compressed symbol data.

    Supports ``key in symb``, ``symb[key]``, and ``get_coeff_from_word``
    without ever materializing the full uncompressed symbol.

    Implementation: each lookup factors into a per-prefix part (heavy, cached
    here) and a per-suffix part (cheap, recomputed). The first query for a
    given prefix computes either ``rest_vals`` (oct, 279 ints from a matrix
    application) or ``indep_values`` (quad, 24 floats); subsequent queries
    for keys with the same prefix skip that and only run the per-suffix
    sparse column dot product against the cached intermediate.

    Critical for the L=7 F/B-anchored sampler workload, where ~1.25M term
    lookups touch on the order of ~1k unique prefixes -- without the cache
    each lookup redoes the matrix work for the prefix; with the cache it's
    paid once per prefix.
    """

    def __init__(self, loop, comp_format="quad", data=None, seeds=None):
        from aiamplitudes_common_public import Phi2Symb

        self.loop = loop
        self.comp_format = comp_format
        if data is not None:
            self.data = data
        else:
            self.data = Phi2Symb(loop, comp_format)
        # Per-prefix cache. Value is a 279-int rest_vals (oct) or a 24-float
        # indep_values numpy array (quad) -- a few KB per prefix at most.
        self._prefix_cache = {}
        # Optional seed-dict fast path: if the key is in here, return that
        # value immediately and skip the prefix-cache + column dot product.
        # Caller-supplied; at L=7/8 this is typically the raydict from
        # extract_seeds (d-rays + a-rays, ~hundreds of thousands of keys).
        self.seeds = seeds if seeds is not None else {}
        if comp_format == "oct":
            self._suffix_len = 8
            _load_oct_matrices()  # pre-load the per-letter B matrices
        else:
            self._suffix_len = 4

    def _coef_oct(self, key):
        prefix = key[:-8]
        suffix = key[-8:]
        m = _oct_matrix_cache[prefix[-1]]
        if suffix not in m['sf_index']:
            return 0
        rv = self._prefix_cache.get(prefix)
        if rv is None:
            rv = _oct_rest_vals(prefix, self.data)
            self._prefix_cache[prefix] = rv
        return _oct_term_from_rest_vals(rv, suffix, prefix[-1])

    def _coef_quad(self, key):
        prefix = key[:-4]
        suffix = key[-4:]
        _, _, sf_index = _get_quad_matrix(prefix[-1])
        if suffix not in sf_index:
            return 0
        iv = self._prefix_cache.get(prefix)
        if iv is None:
            iv = _quad_indep_values(prefix, self.data)
            self._prefix_cache[prefix] = iv
        return _quad_term_from_indep(iv, suffix, prefix[-1])

    def _coef(self, key):
        if len(key) < self._suffix_len:
            return 0
        # Seed-dict fast path: O(1) hit returns immediately, skipping both the
        # prefix-cache and the column dot product entirely.
        v = self.seeds.get(key)
        if v is not None:
            return v
        return self._coef_oct(key) if self.comp_format == "oct" else self._coef_quad(key)

    def __getitem__(self, key):
        return self._coef(key)

    def __contains__(self, key):
        return self._coef(key) != 0

    def get(self, key, default=None):
        v = self._coef(key)
        return v if v != 0 else default

    def __bool__(self):
        return True

    def __repr__(self):
        return (f"LazySymbol(loop={self.loop}, format={self.comp_format!r}, "
                f"cached_prefixes={len(self._prefix_cache)})")

    def clear_cache(self):
        """Drop the per-prefix intermediate cache. Use during a long-running
        iteration if memory grows past comfort."""
        self._prefix_cache.clear()
