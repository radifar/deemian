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
    measurement.interacting_subjects["protein_A:oseltamivir"] = ("protein_A", "oseltamivir")
    measurement.conformation_range("1", "4")
    measurement.conformation.extend([5])

    result = data.calculate_interactions()

    assert (
        repr(measurement)
        == """Measurement(interactions=['all'], \
ionizable={'positive': True, 'negative': False}, \
interacting_subjects={'protein_A:oseltamivir': ('protein_A', 'oseltamivir')}, \
conformation=[1, 2, 3, 4, 5])"""
    )
    assert result == 1


def test_deemian_data_presentation():
    data = DeemianData()

    readable_output = data.write_readable_output("protein_ligand.txt", "protein_ligand")
    deemian_data = data.write_deemian_data("protein_ligand.db", "protein_ligand")

    assert readable_output == ("protein_ligand.txt", "protein_ligand")
    assert deemian_data == ("protein_ligand.db", "protein_ligand")
