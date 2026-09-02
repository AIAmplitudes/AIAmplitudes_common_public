from aiamplitudes_common_public import bspace, fspace
from aiamplitudes_common_public.fbspaces import preload_fbspaces


def test_front_dualspace_is_available_through_fspace():
    dual = fspace(2, "PD")

    assert set(dual) == {f"FPD_2_{i}" for i in range(1, 10)}


def test_back_dualspace_supports_phi2_and_phi3():
    phi2 = bspace(2, "PD", optrace="phi2")
    phi3 = bspace(2, "PD", optrace="phi3")

    assert set(phi2) == {f"BPD_2_{i}" for i in range(1, 7)}
    assert set(phi3) == {f"BPD_2_{i}" for i in range(1, 7)}
    assert phi2 != phi3


def test_preloader_routes_pd_backspace_optrace():
    loaded = preload_fbspaces("back", None, None, "PD", 2, optrace="phi3")

    assert loaded["bspace"] == bspace(2, "PD", optrace="phi3")
