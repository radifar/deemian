from pathlib import Path

from lark import Tree
import pytest

from deemian.engine import processor


@pytest.fixture
def n1_5nzn():
    base_path = Path(__file__).parent

    with open(base_path / "data/n1-oseltamivir-5nzn.txt") as reader:
        text = reader.read()
        return processor.parser(text)


@pytest.fixture
def spike_7u0n():
    base_path = Path(__file__).parent

    with open(base_path / "data/omicron-rbd-ace2-7u0n-multiselect.txt") as reader:
        text = reader.read()
        return processor.parser(text)


@pytest.fixture
def vps4_2kw3():
    base_path = Path(__file__).parent

    with open(base_path / "data/vps4-chmp6-2kw3.txt") as reader:
        text = reader.read()
        return processor.parser(text)


def count_instruction(tree: Tree, data: str):
    instruction_list = list(tree.find_data(data))

    return len(instruction_list)


def test_instruction_count(n1_5nzn, spike_7u0n, vps4_2kw3):
    trees = [n1_5nzn, spike_7u0n, vps4_2kw3]

    select_count = [count_instruction(tree, "assign_selection") for tree in trees]
    bond_correction_count = [count_instruction(tree, "bond_correction") for tree in trees]
    interaction_list_count = [count_instruction(tree, "interaction_list") for tree in trees]
    ionizable_count = [count_instruction(tree, "include_ionizable") for tree in trees]
    interacting_subject_count = [count_instruction(tree, "interacting_subject") for tree in trees]
    interacting_subject_with_alias_count = [
        count_instruction(tree, "interacting_subject_with_alias") for tree in trees
    ]
    conformation_range_count = [count_instruction(tree, "conformation_range") for tree in trees]
    readable_output_count = [count_instruction(tree, "readable_output") for tree in trees]
    deemian_data_count = [count_instruction(tree, "deemian_data") for tree in trees]

    assert select_count == [2, 3, 2]
    assert bond_correction_count == [1, 0, 0]
    assert interaction_list_count == [1, 1, 1]
    assert ionizable_count == [2, 2, 2]
    assert interacting_subject_count == [1, 0, 1]
    assert interacting_subject_with_alias_count == [0, 4, 0]
    assert conformation_range_count == [0, 0, 1]
    assert readable_output_count == [1, 1, 1]
    assert deemian_data_count == [1, 1, 1]
