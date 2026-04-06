"""
Lazy symbol interface for high-loop data.

Provides dict-like access to symbol coefficients without holding the full
uncompressed symbol in memory. Backed by termwise decompression
(UnQuadTerm / UnOctTerm).
"""

from aiamplitudes_common_public.uncompressor import UnQuadTerm, UnOctTerm


class LazySymbol:
    """Dict-like wrapper around compressed symbol data.

    Supports ``key in symb``, ``symb[key]``, and ``get_coeff_from_word``
    without ever materializing the full uncompressed symbol.
    """

    def __init__(self, loop, comp_format="quad", data=None):
        from aiamplitudes_common_public import Phi2Symb

        self.loop = loop
        self.comp_format = comp_format
        if data is not None:
            self.data = data
        else:
            self.data = Phi2Symb(loop, comp_format)
        self._lookup = UnQuadTerm if comp_format == "quad" else UnOctTerm

    def __getitem__(self, key):
        return self._lookup(key, data=self.data)

    def __contains__(self, key):
        return self._lookup(key, data=self.data) != 0

    def get(self, key, default=None):
        v = self._lookup(key, data=self.data)
        return v if v != 0 else default

    def __bool__(self):
        return True

    def __repr__(self):
        return f"LazySymbol(loop={self.loop}, format={self.comp_format!r})"
