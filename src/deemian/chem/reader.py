import pandas as pd


def mol_to_dataframe(mol):
    chain_id = []
    atom_symbol = []
    atom_name = []
    residue_number = []
    residue_name = []
    atom_id = []

    for atom in mol.GetAtoms():
        residue_info = atom.GetPDBResidueInfo()

        atom_id.append(atom.GetIdx())
        atom_symbol.append(atom.GetSymbol())
        atom_name.append(atom.GetMonomerInfo().GetName())
        chain_id.append(residue_info.GetChainId())
        residue_number.append(residue_info.GetResidueNumber())
        residue_name.append(residue_info.GetResidueName())

    # note: string PyArrow is faster, int64 as default for integer is also faster but use more memory
    mol_record = dict(
        chain_id=pd.Series(chain_id, dtype="string[pyarrow]"),
        atom_symbol=pd.Series(atom_symbol, dtype="string[pyarrow]"),
        atom_name=pd.Series(atom_name, dtype="string[pyarrow]"),
        residue_number=pd.Series(residue_number),
        residue_name=pd.Series(residue_name, dtype="string[pyarrow]"),
        atom_id=pd.Series(atom_id),
    )

    # note: using pyarrow for list of float making it slower for 'to_list' conversion,
    #       the scipy's KDTree require list of list as the input.
    for index, conformer in enumerate(mol.GetConformers()):
        coordinate_number = "conf_" + str(index)
        mol_record[coordinate_number] = list(conformer.GetPositions())

    return pd.DataFrame(mol_record).set_index("atom_id")
