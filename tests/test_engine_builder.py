from unittest.mock import Mock, call, patch

from deemian.engine.builder import Chem, DeemianData


def test_deemian_data_molecule_and_selection():
    data = DeemianData()

    data.add_molecule("tests/data/5nzn.pdb")
    data.add_selection("protein_A", [("chain", "A"), ("and", "protein")], "tests/data/5nzn.pdb")
    data.add_selection("oseltamivir", [("chain", "A"), ("and", "resname", "G39")], "tests/data/5nzn.pdb")

    molecule = data.molecule["tests/data/5nzn.pdb"]
    protein_A = data.selection["protein_A"]
    oseltamivir = data.selection["oseltamivir"]

    assert molecule.rdkit_mol.GetNumAtoms() == 27933
    assert molecule.mol_dataframe.shape == (27933, 6)
    assert protein_A.mol_parent == "tests/data/5nzn.pdb"
    assert protein_A.mol_dataframe.shape == (5819, 6)
    assert protein_A.mol_pdb_block is None
    assert oseltamivir.mol_parent == "tests/data/5nzn.pdb"
    assert oseltamivir.mol_dataframe.shape == (44, 6)
    assert oseltamivir.mol_pdb_block is None


def test_deemian_data_correct_bond():
    data = DeemianData()

    data.add_molecule("tests/data/5nzn.pdb")
    data.add_selection("oseltamivir", [("chain", "A"), ("and", "resname", "G39")], "tests/data/5nzn.pdb")
    data.correct_bond("oseltamivir", "CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)O")

    oseltamivir_mol = data.molecule["oseltamivir"]
    oseltamivir_selection = data.selection["oseltamivir"]

    oseltamivir_smi = Chem.MolToSmiles(oseltamivir_mol.rdkit_mol)
    template_smi = "CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)O"

    oseltamivir_canon_smi = Chem.CanonSmiles(oseltamivir_smi)
    template_canon_smi = Chem.CanonSmiles(template_smi)

    assert oseltamivir_canon_smi == template_canon_smi
    assert oseltamivir_selection.mol_parent == "oseltamivir"
    assert oseltamivir_selection.mol_dataframe.shape == (20, 6)
    assert oseltamivir_selection.mol_pdb_block is not None


def test_deemian_data_measurement():
    data = DeemianData()

    measurement = data.add_measurement("protein_ligand")
    measurement.interactions.extend(["all"])
    measurement.set_ionizable("positive", "true")
    measurement.set_ionizable("negative", "false")
    measurement.interacting_subjects["oseltamivir:protein_A"] = ("oseltamivir", "protein_A")
    measurement.set_conformation_range("1", "4")
    measurement.conformation.extend([5])

    with patch.object(data, "calculate_interactions", return_value=1):
        result = data.calculate_interactions("protein_ligand")

    assert (
        repr(measurement)
        == """Measurement(interactions=['all'], \
ionizable={'positive': True, 'negative': False}, \
interacting_subjects={'oseltamivir:protein_A': ('oseltamivir', 'protein_A')}, \
conformation=[1, 2, 3, 4, 5], conformation_range=['1', '4'], \
calculation_results={})"""
    )
    assert result == 1


def test_deemian_data_calculate_interactions():
    data = DeemianData()

    neuraminidase = Mock()
    neuraminidase.rdkit_mol = "neuraminidase:rdkitmol"
    data.molecule["5nzn.pdb"] = neuraminidase

    protein_A = Mock()
    protein_A.mol_parent = "5nzn.pdb"
    protein_A.mol_dataframe = "protein_A:pd.DataFrame"

    oseltamivir = Mock()
    oseltamivir.mol_parent = "5nzn.pdb"
    oseltamivir.mol_dataframe = "oseltamivir:pd.DataFrame"

    data.selection["protein_A"] = protein_A
    data.selection["oseltamivir"] = oseltamivir

    measurement = data.add_measurement("protein_ligand")
    measurement.interactions.extend(["all"])
    measurement.set_ionizable("positive", "true")
    measurement.set_ionizable("negative", "true")
    measurement.interacting_subjects["oseltamivir:protein_A"] = ("oseltamivir", "protein_A")
    measurement.conformation.extend([1])

    with patch("deemian.engine.builder.InteractionData") as int_data:
        data.calculate_interactions("protein_ligand")
        int_data.assert_has_calls(
            [
                call(
                    "neuraminidase:rdkitmol",
                    "neuraminidase:rdkitmol",
                    "oseltamivir:pd.DataFrame",
                    "protein_A:pd.DataFrame",
                    [1],
                ),
                call().calculate_electrostatic(positive=True, negative=True),
            ]
        )


def test_deemian_data_calculate_interactions_empty_conformation():
    data = DeemianData()

    neuraminidase = Mock()
    neuraminidase.rdkit_mol = "neuraminidase:rdkitmol"
    data.molecule["5nzn.pdb"] = neuraminidase

    protein_A = Mock()
    protein_A.mol_parent = "5nzn.pdb"
    protein_A.mol_dataframe = "protein_A:pd.DataFrame"

    oseltamivir = Mock()
    oseltamivir.mol_parent = "5nzn.pdb"
    oseltamivir.mol_dataframe = "oseltamivir:pd.DataFrame"

    data.selection["protein_A"] = protein_A
    data.selection["oseltamivir"] = oseltamivir

    measurement = data.add_measurement("protein_ligand")
    measurement.interactions.extend(["all"])
    measurement.set_ionizable("positive", "true")
    measurement.set_ionizable("negative", "true")
    measurement.interacting_subjects["oseltamivir:protein_A"] = ("oseltamivir", "protein_A")

    with patch("deemian.engine.builder.InteractionData") as int_data:
        data.calculate_interactions("protein_ligand")
        int_data.assert_called_once()


def test_deemian_data_presentation():
    data = DeemianData()

    protein_A = Mock()
    protein_A.mol_parent = "5nzn.pdb"
    protein_A.mol_dataframe = "protein_A:pd.DataFrame"
    protein_A.mol_pdb_block = None

    oseltamivir = Mock()
    oseltamivir.mol_parent = "oseltamivir_corrected.pdb"
    oseltamivir.mol_dataframe = "oseltamivir:pd.DataFrame"
    oseltamivir.mol_pdb_block = "pdb_block:str"

    data.selection["protein_A"] = protein_A
    data.selection["oseltamivir"] = oseltamivir

    measurement = data.add_measurement("protein_ligand")
    measurement.interactions.extend(["all"])
    measurement.set_ionizable("positive", "true")
    measurement.set_ionizable("negative", "true")
    measurement.interacting_subjects["oseltamivir:protein_A"] = ("oseltamivir", "protein_A")

    interaction_data = Mock()
    interaction_data.dataframe = "interaction_data:pd.DataFrame"
    measurement.calculation_results["oseltamivir:protein_A"] = interaction_data

    with patch("deemian.engine.builder.generate_report") as reporter:
        with patch("deemian.engine.builder.write_readable") as writer:
            data.write_readable_output("protein_ligand", "protein_ligand.txt", "detailed_conf_first")

            reporter.assert_called_once()
            writer.assert_called_once()

    with patch("deemian.engine.builder.write_corrected_molecule") as writer_1:
        with patch("deemian.engine.builder.write_calculation_result") as writer_2:
            with patch("deemian.engine.builder.write_metadata") as metadata_writer:
                with patch("deemian.engine.builder.write_bundle") as bundler:
                    data.write_deemian_data("protein_ligand", "protein_ligand.dd")

                    writer_1.assert_called_once()
                    writer_2.assert_called_once()
                    metadata_writer.assert_called_once()
                    bundler.assert_called_once()
