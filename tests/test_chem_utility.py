import pandas as pd
from rdkit.Chem import AllChem as Chem

from deemian.chem.utility import dataframe_to_pdb_block


def test_dataframe_to_pdb_block():
    """
    Test dataframe_to_pdb_block using oseltamivir dataframe, which is extracted
    from 5nzn.pdb that loaded with deemian.chem.reader.mol_to_dataframe and
    selected where chain == A and resname == G39, then saved into parquet file
    """
    oseltamivir_df = pd.read_parquet("tests/data/oseltamivir.parquet.gzip")

    oseltamivir_pdb = dataframe_to_pdb_block(oseltamivir_df)
    oseltamivir_mol = Chem.MolFromPDBBlock(oseltamivir_pdb)

    oseltamivir_atom_num = oseltamivir_mol.GetNumAtoms()

    assert oseltamivir_atom_num == 20
