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


@pytest.fixture
def n1_transformed(n1_5nzn):
    instruction_transformer = processor.InstructionTransformer()

    return instruction_transformer.transform(n1_5nzn)


@pytest.fixture
def spike_transformed(spike_7u0n):
    instruction_transformer = processor.InstructionTransformer()

    return instruction_transformer.transform(spike_7u0n)


@pytest.fixture
def vps_transformed(vps4_2kw3):
    instruction_transformer = processor.InstructionTransformer()

    return instruction_transformer.transform(vps4_2kw3)


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


def test_molecule_transformer(n1_transformed, spike_transformed, vps_transformed):
    n1_molecule = n1_transformed.children[0]
    n1_selection_1 = n1_molecule.children[1]
    n1_selection_2 = n1_molecule.children[2]
    n1_bond_correction = n1_molecule.children[3]

    spike_molecule = spike_transformed.children[0]
    spike_selection_1 = spike_molecule.children[1]
    spike_selection_2 = spike_molecule.children[2]
    spike_selection_3 = spike_molecule.children[3]

    vps_molecule = vps_transformed.children[0]
    vps_selection_1 = vps_molecule.children[1]
    vps_selection_2 = vps_molecule.children[2]

    assert n1_selection_1.name == "protein_A"
    assert n1_selection_1.selection == [("chain", "A"), ("and", "protein")]
    assert n1_selection_1.type == "selection"
    assert n1_selection_2.name == "oseltamivir"
    assert n1_selection_2.selection == [("chain", "A"), ("and", "resname", "G39")]
    assert n1_selection_2.type == "selection"

    assert spike_selection_1.name == "ace2_hotspot31"
    assert spike_selection_1.selection == [("chain", "A"), ("and", "protein"), ("and", "resid_range", "1", "55")]
    assert spike_selection_1.type == "selection"
    assert spike_selection_2.name == "ace2_hotspot353"
    assert spike_selection_2.selection == [("chain", "A"), ("and", "protein"), ("and", "resid_range", "341", "364")]
    assert spike_selection_2.type == "selection"
    assert spike_selection_3.name == "spike_rbm"
    assert spike_selection_3.selection == [("chain", "E"), ("and", "protein"), ("and", "resid_range", "470", "510")]
    assert spike_selection_3.type == "selection"

    assert vps_selection_1.name == "vps4"
    assert vps_selection_1.selection == [("protein",), ("and", "chain", "A")]
    assert vps_selection_1.type == "selection"
    assert vps_selection_2.name == "chmp6"
    assert vps_selection_2.selection == [("protein",), ("and", "chain", "B")]
    assert vps_selection_2.type == "selection"

    assert n1_bond_correction.name == "oseltamivir"
    assert n1_bond_correction.template == "CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)O"
    assert n1_bond_correction.type == "bond_correction"


def test_measure_transform(n1_transformed, spike_transformed, vps_transformed):
    n1_measure = n1_transformed.children[1]
    n1_interactions = n1_measure.children[1]
    n1_ionizable1 = n1_measure.children[2]
    n1_ionizable2 = n1_measure.children[3]
    n1_interacting_subject = n1_measure.children[4]
    n1_conformation = n1_measure.children[5]

    spike_measure = spike_transformed.children[1]
    spike_interacting_subject_as = spike_measure.children[4]

    vps_measure = vps_transformed.children[1]
    vps_conformation_range = vps_measure.children[5]

    assert n1_interactions.interactions == ["all"]
    assert n1_interactions.type == "interaction_list"
    assert n1_ionizable1.charge == "positive"
    assert n1_ionizable1.boolean == "true"
    assert n1_ionizable1.type == "include_ionizable"
    assert n1_ionizable2.charge == "negative"
    assert n1_ionizable2.boolean == "true"
    assert n1_ionizable2.type == "include_ionizable"
    assert n1_interacting_subject.subject_1 == "protein_A"
    assert n1_interacting_subject.subject_2 == "oseltamivir"
    assert n1_interacting_subject.name == "protein_A:oseltamivir"
    assert n1_interacting_subject.type == "interacting_subject"
    assert n1_conformation.number == "1"
    assert n1_conformation.type == "conformation"

    assert spike_interacting_subject_as.subject_1 == "spike_rbm"
    assert spike_interacting_subject_as.subject_2 == "ace2_hotspot31"
    assert spike_interacting_subject_as.name == "rbm_h31"
    assert spike_interacting_subject_as.type == "interacting_subject"

    assert vps_conformation_range.start == "1"
    assert vps_conformation_range.end == "20"
    assert vps_conformation_range.type == "conformation_range"


def test_present_transform(n1_transformed):
    n1_present = n1_transformed.children[2]
    n1_readable_output = n1_present.children[1]
    n1_deemian_data = n1_present.children[2]

    assert n1_readable_output.results == ["interactions"]
    assert n1_readable_output.out_file == "protein_ligand.txt"
    assert n1_readable_output.type == "readable_output"
    assert n1_deemian_data.out_file == "protein_ligand.db"
    assert n1_deemian_data.type == "deemian_data"
