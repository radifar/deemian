from rdkit.Chem import AllChem as Chem

from deemian.chem.utility import dataframe_to_pdb_block, flatten_atom_list


def test_flatten_atom_list():
    test_list = [[1, 2, 3], [4, 5, 6], [1, 10, 12]]
    expected_list = [1, 2, 3, 4, 5, 6, 10, 12]

    result_list = flatten_atom_list(test_list)

    assert result_list == expected_list


def test_dataframe_to_pdb_block(oseltamivir_df):
    """
    Test dataframe_to_pdb_block using oseltamivir dataframe, which is extracted
    from 5nzn.pdb that loaded with deemian.chem.reader.mol_to_dataframe and
    selected where chain == A and resname == G39, then saved into parquet file
    """

    oseltamivir_pdb = dataframe_to_pdb_block(oseltamivir_df)
    oseltamivir_mol = Chem.MolFromPDBBlock(oseltamivir_pdb)

    oseltamivir_atom_num = oseltamivir_mol.GetNumAtoms()

    assert oseltamivir_atom_num == 20
