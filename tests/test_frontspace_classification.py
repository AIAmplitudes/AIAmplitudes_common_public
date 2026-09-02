from aiamplitudes_common_public import (
    br_dep,
    br_indep,
    br_mand,
    br_rels,
    bspace,
    fr_dep,
    fr_indep,
    fr_mand,
    fr_rels,
    fspace,
)
from aiamplitudes_common_public.fbspaces import F_number, NF_rels


def test_frontspace_classification_at_weight_two():
    assert fr_mand(2) == {
        "FR_2_1": "aa",
        "FR_2_2": "bb",
        "FR_2_3": "cc",
    }
    assert fr_indep(2) == {
        "FR_2_4": "bc",
        "FR_2_5": "ca",
        "FR_2_6": "ab",
        "FR_2_7": "bd",
        "FR_2_8": "ce",
        "FR_2_9": "af",
    }
    assert fr_dep(2) == fr_rels(2)


def test_mandatory_and_independent_partition_restrictive_basis():
    for weight in range(1, 5):
        mandatory = fr_mand(weight)
        independent = fr_indep(weight)

        assert mandatory.keys().isdisjoint(independent)
        assert mandatory | independent == fspace(weight, "R")


def test_dependents_are_disjoint_from_restrictive_basis_words():
    for weight in range(1, 5):
        assert set(fr_dep(weight)).isdisjoint(fspace(weight, "R").values())


def test_mandatory_and_independent_partition_restrictive_back_basis():
    for weight in range(1, 5):
        mandatory = br_mand(weight)
        independent = br_indep(weight)

        assert mandatory.keys().isdisjoint(independent)
        assert mandatory | independent == bspace(weight, "R")


def test_back_dependents_match_relations_and_are_disjoint_from_basis_words():
    for weight in range(1, 5):
        assert br_dep(weight) == br_rels(weight)
        assert set(br_dep(weight)).isdisjoint(bspace(weight, "R").values())


def test_front_dimension_tables_match_available_data():
    for weight in range(1, 8):
        assert F_number[weight] == len(fspace(weight, "R"))
        assert NF_rels[weight] == len(fr_rels(weight))

    assert F_number[8] is None
    assert NF_rels[8] is None
