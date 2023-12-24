from itertools import chain

import pandas as pd


def flatten_atom_list(nested_list_of_atoms):
    return list(set(chain.from_iterable(nested_list_of_atoms)))


def dataframe_to_pdb_block(df: pd.DataFrame) -> str:
    df_noh = df[df["atom_symbol"] != "H"]
    # from https://stackoverflow.com/questions/35491274/split-a-pandas-column-of-lists-into-multiple-columns

    df_pdb = pd.DataFrame(df_noh["conf_1"].to_list(), columns=["x", "y", "z"])
    df_pdb.index = df_noh.index
    df_pdb = df_pdb.join(df_noh).reset_index()
    df_pdb["record_type"] = "HETATM"
    df_pdb["occupancy"] = "1.00"
    df_pdb["temperature"] = "0.00"
    df_pdb = df_pdb[
        [
            "record_type",
            "atom_id",
            "atom_name",
            "residue_name",
            "chain_id",
            "residue_number",
            "x",
            "y",
            "z",
            "occupancy",
            "temperature",
            "atom_symbol",
        ]
    ]
    pdb_string = df_pdb.to_string(index=False, header=False)

    # https://stackoverflow.com/questions/17796017/how-do-i-output-a-pdb-file-using-python-script
    pdb_lines = []
    for line in pdb_string.splitlines():
        values = line.split()
        values[0] = values[0].ljust(6)  # atom#6s
        values[1] = values[1].rjust(5)  # atomnum#5d
        values[2] = values[2].center(4)  # atomname#4s
        values[3] = values[3].ljust(3)  # resname#1s
        values[4] = values[4].rjust(1)  # Astring
        values[5] = values[5].rjust(4)  # resnum
        values[6] = str("%8.3f" % (float(values[6]))).rjust(8)  # x
        values[7] = str("%8.3f" % (float(values[7]))).rjust(8)  # y
        values[8] = str("%8.3f" % (float(values[8]))).rjust(8)  # z
        values[9] = str("%6.2f" % (float(values[9]))).rjust(6)  # occ
        values[10] = str("%6.2f" % (float(values[10]))).ljust(6)  # temp
        values[11] = values[11].rjust(12)  # elname
        pdb_lines.append("%s%s %s %s %s%s    %s%s%s%s%s%s" % tuple(values))
    pdb_block = "\n".join(pdb_lines)

    return pdb_block
