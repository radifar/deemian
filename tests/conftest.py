import pytest
import pandas as pd

from rdkit import Chem

from deemian.chem.selection import mol_dataframe_selection


@pytest.fixture
def n1():
    return Chem.MolFromPDBFile("tests/data/5nzn.pdb", removeHs=False)


@pytest.fixture
def oseltamivir_corrected():
    return Chem.MolFromPDBFile("tests/data/oseltamivir_corrected.pdb", removeHs=False)


@pytest.fixture
def n1_df():
    """
    Parquet file originated from 5nzn.pdb, protonated using Protoss.
    Loaded with deemian.chem.reader.mol_to_dataframe saved using pyarrow engine
    """

    return pd.read_parquet("tests/data/5nzn.parquet.gzip")


@pytest.fixture
def spike_df():
    """
    Parquet file originated from 7u0n.pdb, protonated using Protoss.
    Loaded with deemian.chem.reader.mol_to_dataframe saved using pyarrow engine
    """

    return pd.read_parquet("tests/data/5nzn.parquet.gzip")


@pytest.fixture
def vps4_df():
    """
    Parquet file originated from 2k3w.pdb without preparation (it is NMR solution structure)
    Loaded with deemian.chem.reader.mol_to_dataframe saved using pyarrow engine
    """

    return pd.read_parquet("tests/data/2k3w.parquet.gzip")


@pytest.fixture
def n1_protein_A(n1_df):
    protein_A_selection = [("chain", "A"), ("and", "protein")]

    return mol_dataframe_selection(protein_A_selection, n1_df)


@pytest.fixture
def oseltamivir_df():
    return pd.read_parquet("tests/data/oseltamivir.parquet.gzip")


@pytest.fixture
def oseltamivir_corrected_df():
    return pd.read_parquet("tests/data/oseltamivir_corrected.parquet.gzip")
