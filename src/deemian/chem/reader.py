import pandas as pd


def mol_to_dataframe(mol):
    mol_record = dict(chain_id=[], atom_symbol=[], atom_name=[], residue_number=[], residue_name=[], atom_id=[])

    for atom in mol.GetAtoms():
        residue_info = atom.GetPDBResidueInfo()

        atom_id = atom.GetIdx()
        atom_symbol = atom.GetSymbol()
        atom_name = atom.GetMonomerInfo().GetName()
        chain_id = residue_info.GetChainId()
        residue_number = residue_info.GetResidueNumber()
        residue_name = residue_info.GetResidueName()

        mol_record["atom_id"].append(atom_id)
        mol_record["atom_symbol"].append(atom_symbol)
        mol_record["atom_name"].append(atom_name)
        mol_record["chain_id"].append(chain_id)
        mol_record["residue_number"].append(residue_number)
        mol_record["residue_name"].append(residue_name)

    for index, conformer in enumerate(mol.GetConformers()):
        coordinate_number = "conf_" + str(index)
        mol_record[coordinate_number] = list(conformer.GetPositions())

    return pd.DataFrame(mol_record).set_index("atom_id")
