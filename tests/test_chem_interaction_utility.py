import pytest
from rdkit.Chem import AllChem as Chem

from deemian.chem.interaction_utility import filter_charge, generate_pair_info


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


def test_generate_pair_info(n1, oseltamivir_corrected, n1_protein_A, oseltamivir_corrected_df):
    g39_all_cation = filter_charge(oseltamivir_corrected, "all_cation")
    g39_all_anion = filter_charge(oseltamivir_corrected, "all_anion")
    n1_all_cation = filter_charge(n1, "all_cation")
    n1_all_anion = filter_charge(n1, "all_anion")

    exclude_atom = ["C", "N", "S", "P"]
    cation_1_df = oseltamivir_corrected_df[oseltamivir_corrected_df.index.isin(g39_all_cation)]
    anion_1_df = oseltamivir_corrected_df[
        oseltamivir_corrected_df.index.isin(g39_all_anion)
        & ~oseltamivir_corrected_df["atom_symbol"].isin(exclude_atom)
    ]
    cation_2_df = n1_protein_A[n1_protein_A.index.isin(n1_all_cation)]
    anion_2_df = n1_protein_A[n1_protein_A.index.isin(n1_all_anion) & ~n1_protein_A["atom_symbol"].isin(exclude_atom)]

    anion_nearby = [[4, 5, 10]]
    cation_nearby = [[9, 10, 11, 84, 85, 86], [59, 60, 61, 84, 85, 86]]
    pair_info_as_cation = generate_pair_info(anion_nearby, cation_1_df, anion_2_df, 1, "electrostatic")
    pair_info_as_anion = generate_pair_info(cation_nearby, anion_1_df, cation_2_df, 1, "electrostatic")

    assert pair_info_as_cation.shape == (3, 17)
    assert pair_info_as_anion.shape == (12, 17)
    assert bool((pair_info_as_cation["distance"] < 4.5).all()) is True
    assert bool((pair_info_as_anion["distance"] < 4.5).all()) is True
