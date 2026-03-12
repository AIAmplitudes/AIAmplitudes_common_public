"""
Tests for uncompressor.py — verifies UnQuad and UnQuadLoop produce results
matching the full (uncompressed) symbol from Phi2Symb.
"""
from aiamplitudes_common_public import Phi2Symb, UnQuad, UnQuadLoop


def compare(unq, corr, label=""):
    unq_only = [(k, v) for k, v in unq.items() if k not in corr]
    diff = [(k, unq[k], corr[k]) for k in unq if k in corr and unq[k] != corr[k]]
    corr_only = [(k, v) for k, v in corr.items() if k not in unq]
    print(f"  {label}: unq_only={len(unq_only)}, diff={len(diff)}, corr_only={len(corr_only)}")
    if unq_only:
        print(f"    unq_only samples: {unq_only[:5]}")
    if diff:
        print(f"    diff samples: {diff[:5]}")
    if corr_only:
        print(f"    corr_only samples: {corr_only[:5]}")
    return len(unq_only) == 0 and len(diff) == 0 and len(corr_only) == 0


def test_unquadloop(loops=None):
    if loops is None:
        loops = [3, 4, 5, 6]
    all_pass = True
    for L in loops:
        print(f"=== Loop {L} ===")
        unq = UnQuadLoop(L)
        corr = Phi2Symb(L)
        ok = compare(unq, corr, f"Loop {L}")
        all_pass = all_pass and ok
    return all_pass


def test_unquad_single_prefix():
    # At loop 5, prefix length is 2*5-4 = 6 characters
    prefix = "aabbcd"
    print(f"=== Single prefix: UnQuad('{prefix}') at loop 5 ===")
    data = Phi2Symb(5, "quad")
    unq = UnQuad(prefix, data=data)
    corr = {k: v for k, v in Phi2Symb(5).items() if k.startswith(prefix)}
    return compare(unq, corr, f"{prefix}@L5")


if __name__ == "__main__":
    ok1 = test_unquadloop()
    ok2 = test_unquad_single_prefix()
    print()
    if ok1 and ok2:
        print("All tests passed.")
    else:
        print("SOME TESTS FAILED.")
        exit(1)
