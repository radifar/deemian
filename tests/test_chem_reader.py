from rdkit.Chem import AllChem as Chem

from deemian.chem.reader import mol_to_dataframe


def test_mol_to_dataframe():
    mol = Chem.MolFromPDBFile("tests/data/5nzn.pdb", removeHs=False)
    mol_df = mol_to_dataframe(mol)
    df_column = list(mol_df.columns)

    assert mol_df.shape == (27933, 6)
    assert mol_df.index.name == "atom_id"
    assert df_column == ["chain_id", "atom_symbol", "atom_name", "residue_number", "residue_name", "conf_0"]
