import pandas as pd


def mol_dataframe_selection(selections: list[tuple], mol_df: pd.DataFrame):
    amino_acids = [  # noqa: F841
        "ALA",
        "CYS",
        "ASP",
        "GLU",
        "PHE",
        "GLY",
        "HIS",
        "ILE",
        "LYS",
        "LEU",
        "MET",
        "ASN",
        "PRO",
        "GLN",
        "ARG",
        "SER",
        "THR",
        "VAL",
        "TRP",
        "TYR",
        "PYL",
        "SEC",
        "NVA",
    ]

    query_string = []
    for selection in selections:
        selection = list(selection)
        if selection[0] == "and":
            query_string.append(selection.pop(0))
        if selection[0] == "not":
            query_string.append(selection.pop(0))
        if selection[0] == "protein":
            query_string.append("mol_df.residue_name.isin(amino_acids)")
        elif selection[0] == "resname":
            query_string.append(f"mol_df.residue_name == '{selection[1]}'")
        elif selection[0] == "chain":
            query_string.append(f"mol_df.chain_id == '{selection[1]}'")
        elif selection[0] == "resid":
            ids = [int(id) for id in selection[1:]]  # noqa: F841
            query_string.append("mol_df.residue_number.isin(ids)")
        elif selection[0] == "resid_range":
            start = int(selection[1])
            end = int(selection[2])
            id_range = range(start, end + 1)  # noqa: F841
            query_string.append("mol_df.residue_number.isin(id_range)")

    query_string = " ".join(query_string)

    return mol_df[pd.eval(query_string)]
