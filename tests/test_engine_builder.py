from unittest.mock import patch

import pytest

from deemian.engine.builder import Chem, DeemianData, DeemianDataBuilder


@pytest.fixture
def builder():
    deemian_data = DeemianData()
    deemian_data_builder = DeemianDataBuilder(deemian_data)

    return deemian_data_builder


def test_deemian_data():
    deemian_data = DeemianData()

    assert deemian_data.molecule_dataframe == {}
    assert deemian_data.selections == {}
    assert deemian_data.interactions == []
    assert deemian_data.ionizable == {"positive": False, "negative": False}
    assert deemian_data.interacting_subjects == {}
    assert deemian_data.conformation == []
    assert deemian_data.interaction_details == {}
    assert deemian_data.readable_output == {}


def test_deemian_data_builder_molecule_reader(builder):
    with patch.object(Chem, "MolFromPDBFile") as mol_from_pdb:
        with patch("deemian.engine.builder.mol_to_dataframe") as mol_to_df:
            mol_from_pdb.return_value = "rdkit.Chem.rdchem.Mol"
            mol_to_df.return_value = "pd.DataFrame"

            mol_filename = "tests/data/5nzn.pdb"
            builder.read_molecule(mol_filename)

            deemian_data = builder.generate_deemian_data()

            mol_from_pdb.assert_called_once_with(mol_filename, removeHs=False)
            mol_to_df.assert_called_once_with("rdkit.Chem.rdchem.Mol")
            assert deemian_data.molecule_dataframe[mol_filename] == "pd.DataFrame"


def test_deemian_data_builder_molecule_selection(builder):
    with patch("deemian.engine.builder.mol_dataframe_selection") as mds:
        mds.return_value = "selected_df:pd.DataFrame"
        deemian_data = builder.generate_deemian_data()
        deemian_data.molecule_dataframe["5nzn.pdb"] = "n1_df:pd.DataFrame"

        name = "protein_A"
        selection = [("chain", "A"), ("and", "protein")]
        mol_filename = "5nzn.pdb"

        builder.assign_selection(name, selection, mol_filename)
        correct_bond = builder.correct_bond("oseltamivir", "CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)O")

        mds.assert_called_once_with(selection, "n1_df:pd.DataFrame")
        assert deemian_data.molecule_dataframe["protein_A"] == "selected_df:pd.DataFrame"
        assert deemian_data.selections["protein_A"] == ["5nzn.pdb", ("chain", "A"), ("and", "protein")]
        assert correct_bond == ("oseltamivir", "CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)O")


def test_deemian_data_builder_measure(builder):
    builder.set_interactions(["all"])
    builder.set_ionizable("positive", "true")
    builder.set_ionizable("negative", "false")
    builder.set_interacting_subjects("protein_A", "oseltamivir", "protein_A:oseltamivir")
    builder.set_conformation("1")

    interaction_results = builder.calculate_interactions("protein_ligand")

    deemian_data = builder.generate_deemian_data()

    assert deemian_data.interactions == ["all"]
    assert deemian_data.ionizable["positive"] is True
    assert deemian_data.ionizable["negative"] is False
    assert deemian_data.interacting_subjects["protein_A:oseltamivir"] == ("protein_A", "oseltamivir")
    assert deemian_data.conformation == [1]

    assert interaction_results == "protein_ligand"


def test_deemian_data_builder_measure_conf_range(builder):
    builder.set_conformation_range("1", "10")

    deemian_data = builder.generate_deemian_data()

    assert deemian_data.conformation == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]


def test_deemian_data_builder_present(builder):
    readable_output = builder.write_readable_output("protein_ligand.txt", "protein_ligand")
    deemian_output = builder.write_deemian_data("protein_ligand.db", "protein_ligand")

    assert readable_output == ("protein_ligand.txt", "protein_ligand")
    assert deemian_output == ("protein_ligand.db", "protein_ligand")
