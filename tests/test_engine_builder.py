import pytest

from deemian.engine.builder import DeemianData, DeemianDataBuilder


@pytest.fixture
def builder():
    deemian_data = DeemianData()
    deemian_data_builder = DeemianDataBuilder(deemian_data)

    return deemian_data_builder


def test_deemian_data():
    deemian_data = DeemianData()

    assert deemian_data.selections == {}
    assert deemian_data.interactions == []
    assert deemian_data.ionizable == {"positive": False, "negative": False}
    assert deemian_data.interacting_subjects == {}
    assert deemian_data.conformation == []
    assert deemian_data.interaction_details == {}
    assert deemian_data.readable_output == {}


def test_deemian_data_builder_molecule(builder):
    name = "protein_A"
    selection = [("chain", "A"), ("and", "protein")]
    molecule = "5nzn.pdb"
    builder.assign_selection(name, selection, molecule)
    correct_bond = builder.correct_bond("oseltamivir", "CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)O")

    deemian_data = builder.generate_deemian_data()

    assert deemian_data.selections["protein_A"] == ["5nzn.pdb", ("chain", "A"), ("and", "protein")]
    assert correct_bond == ("oseltamivir", "CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)O")


def test_deemian_data_builder_measure(builder):
    builder.set_interactions(["all"])
    builder.set_ionizable("positive", "true")
    builder.set_ionizable("negative", "false")
    builder.set_interacting_subjects("protein_A", "oseltamivir", "protein_A:oseltamivir")
    builder.set_conformation("1")

    deemian_data = builder.generate_deemian_data()

    assert deemian_data.interactions == ["all"]
    assert deemian_data.ionizable["positive"] is True
    assert deemian_data.ionizable["negative"] is False
    assert deemian_data.interacting_subjects["protein_A:oseltamivir"] == ("protein_A", "oseltamivir")
    assert deemian_data.conformation == [1]


def test_deemian_data_builder_measure_conf_range(builder):
    builder.set_conformation_range("1", "10")

    deemian_data = builder.generate_deemian_data()

    assert deemian_data.conformation == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]


def test_deemian_data_builder_present(builder):
    interaction_results = builder.calculate_interactions()
    readable_output = builder.generate_readable_output()

    assert interaction_results == 1
    assert readable_output == 1
