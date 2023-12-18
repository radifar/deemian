from rdkit.Chem import AllChem as Chem

from deemian.chem.utility import flatten_atom_list


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
