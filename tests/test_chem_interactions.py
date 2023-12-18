import pytest
from rdkit.Chem import AllChem as Chem

from deemian.chem.interactions import filter_charge


@pytest.mark.parametrize(
    "smiles, n_cation, n_anion, n_ionizable_cation, n_ionizable_anion, all_cation, all_anion",
    [
        ("CC(=O)NCCCS(O)(=O)=O", 0, 0, 0, 4, 0, 4),  # acamprosate
        ("CC(=O)NCCCS(=O)(=O)[O-]", 0, 4, 0, 4, 0, 4),  # acamprosate -
        ("CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)O", 0, 0, 1, 3, 1, 3),  # oseltamivir
        ("CCC(CC)O[C@@H]1C=C(C(=O)O)C[C@H]([NH3+])[C@H]1NC(C)=O", 1, 0, 0, 3, 1, 3),  # oseltamivir +
        ("CCC(CC)O[C@@H]1C=C(C(=O)[O-])C[C@H](N)[C@H]1NC(C)=O", 0, 3, 1, 3, 1, 3),  # oseltamivir -
        ("CCC(CC)O[C@@H]1C=C(C(=O)[O-])C[C@H]([NH3+])[C@H]1NC(C)=O", 1, 3, 0, 3, 1, 3),  # oseltamivir +-
        ("NC[C@@H](F)CP(O)=O", 0, 0, 1, 3, 1, 3),  # lesogaberan
        ("[NH3+]C[C@@H](F)C[PH](=O)O", 1, 0, 0, 3, 1, 3),  # lesogaberan +
        ("NC[C@@H](F)C[PH](=O)[O-]", 0, 3, 1, 3, 1, 3),  # lesogaberan -
        ("[NH3+]C[C@@H](F)C[PH](=O)[O-]", 1, 3, 0, 3, 1, 3),  # lesogaberan +-
        ("CN(C)C(=N)NC(N)=N", 0, 0, 7, 0, 7, 0),  # metformin
        ("CN(C)C(=[NH2+])N=C(N)N", 3, 0, 7, 0, 7, 0),  # metformin +
        ("N[C@@H](CCCNC(=N)N[N+]([O-])=O)C(O)=O", 0, 3, 5, 6, 5, 6),  # nitroarginine
        ("NC(=N[N+](=O)[O-])NCCC[C@H]([NH3+])C(=O)O", 1, 3, 4, 6, 5, 6),  # nitroarginine +
        ("NC(=N[N+](=O)[O-])NCCC[C@H](N)C(=O)[O-]", 0, 6, 5, 6, 5, 6),  # nitroarginine -
        ("NC(=N[N+](=O)[O-])NCCC[C@H]([NH3+])C(=O)[O-]", 1, 6, 4, 6, 5, 6),  # nitroarginine +-
        (
            "N=C(N)NCCC[C@H](NC(=O)[C@@H](N)CCCCN)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](Cc1c[nH]cn1)C(=O)O",  # noqa: E501
            0,
            0,
            11,
            9,
            11,
            9,
        ),  # pentapeptide KRDEH generated with RDKit
    ],
)
def test_filter_charge(smiles, n_cation, n_anion, n_ionizable_cation, n_ionizable_anion, all_cation, all_anion):
    """
    test SMILES retrieved from Drugbank, and the charged version generated using SwissParam
    with MMFF+Match parameter. Except Oseltamivir SMILES which retrieved from PDB with ID: G39.
    """
    mol = Chem.MolFromSmiles(smiles)
    n_cation_result = len(filter_charge(mol, "apparent_cation"))
    n_anion_result = len(filter_charge(mol, "apparent_anion"))
    n_ionizable_cation_result = len(filter_charge(mol, "potential_cation"))
    n_ionizable_anion_result = len(filter_charge(mol, "potential_anion"))
    all_cation_result = len(filter_charge(mol, "all_cation"))
    all_anion_result = len(filter_charge(mol, "all_anion"))

    assert n_cation == n_cation_result
    assert n_anion == n_anion_result
    assert n_ionizable_cation == n_ionizable_cation_result
    assert n_ionizable_anion == n_ionizable_anion_result
    assert all_cation == all_cation_result
    assert all_anion == all_anion_result


def test_filter_charge_keyerror():
    with pytest.raises(KeyError) as e_info:
        mol = Chem.MolFromSmiles("CC(=O)NCCCS(O)(=O)=O")
        filter_charge(mol, "cation")
    assert str(e_info.value) == "'Error: filtering mode cation is not recognized'"
