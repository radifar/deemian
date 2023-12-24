from numpy.linalg import norm as dist
import pandas as pd
from rdkit.Chem import AllChem as Chem

from deemian.chem.utility import flatten_atom_list


def generate_pair_info(query_result, s1_df, s2_df, conformation, interaction_type):
    pair_info_df = pd.DataFrame()
    for i, s2_list in enumerate(query_result):
        s1 = s1_df.iloc[[i]].reset_index()
        s2 = s2_df.iloc[s2_list].reset_index()
        s2["atom_id_1"] = s1.iloc[0].atom_id
        pair = s1.merge(s2, left_on="atom_id", right_on="atom_id_1", suffixes=("_s1", "_s2")).drop(
            columns=["atom_id_1"]
        )
        pair_info_df = pd.concat([pair_info_df, pair])

    pair_info_df.reset_index(drop=True, inplace=True)

    conf_1 = f"conf_{conformation}_s1"
    conf_2 = f"conf_{conformation}_s2"
    pair_info_df["distance"] = (pair_info_df[conf_1] - pair_info_df[conf_2]).map(dist)
    pair_info_df["conformation"] = conformation
    pair_info_df["interaction_type"] = interaction_type

    return pair_info_df


def filter_charge(mol: Chem.rdchem.Mol, mode: str) -> list[int]:
    apparent_cation = ["[$([*+1,*+2,*+3]);!$([N+]-[O-])]", "[NX3&!$([NX3]-O)]-[C]=[NX3+]"]
    apparent_anion = ["O=[C,S,N,P]-[O-]", "[*-1,*-2]"]
    apparent_cation = [Chem.MolFromSmarts(p) for p in apparent_cation]  # p for pattern
    apparent_anion = [Chem.MolFromSmarts(p) for p in apparent_anion]

    potential_cation = [
        "[$([N;H2&+0][C;!$(C=*)]);!$(N[a])]",
        "[$([N;H1&+0]([C;!$(C=*)])[C;!$(C=*)]);!$(N[a])]",
        "[$([N;H0&+0]([C;!$(C=*)])([C;!$(C=*)])[C;!$(C=*)]);!$(N[a])]",
        "NC(=N)",
        "[n;R1]1[c;R1][n;R1][c;R1][c;R1]1",
    ]
    potential_anion = ["O=[C,S,N,P]-[OH,O-]"]
    potential_cation = [Chem.MolFromSmarts(p) for p in potential_cation]
    potential_anion = [Chem.MolFromSmarts(p) for p in potential_anion]

    all_cation = apparent_cation + potential_cation
    all_anion = apparent_anion + potential_anion

    filter_dictionary = dict(
        apparent_cation=apparent_cation,
        apparent_anion=apparent_anion,
        potential_cation=potential_cation,
        potential_anion=potential_anion,
        all_cation=all_cation,
        all_anion=all_anion,
    )

    combined_match = []
    try:
        for pattern in filter_dictionary[mode]:
            matched = flatten_atom_list(mol.GetSubstructMatches(pattern))
            combined_match.append(matched)
    except KeyError:
        raise KeyError(f"Error: filtering mode {mode} is not recognized")

    return flatten_atom_list(combined_match)
