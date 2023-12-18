from dataclasses import dataclass, field

from rdkit.Chem import AllChem as Chem

from deemian.chem.reader import mol_to_dataframe
from deemian.chem.selection import mol_dataframe_selection
from deemian.chem.utility import dataframe_to_pdb_block


@dataclass
class DeemianData:
    molecule_dataframe: dict = field(default_factory=lambda: {})
    molecule_pdb_block: dict = field(default_factory=lambda: {})
    selections: dict = field(default_factory=lambda: {})
    interactions: list = field(default_factory=lambda: [])
    ionizable: dict = field(default_factory=lambda: {"positive": False, "negative": False})
    interacting_subjects: dict = field(default_factory=lambda: {})
    conformation: list = field(default_factory=lambda: [])
    interaction_details: dict = field(default_factory=lambda: {})
    readable_output: dict = field(default_factory=lambda: {})


class DeemianDataBuilder:
    def __init__(self, deemian_data: DeemianData) -> None:
        self.deemian_data = deemian_data

    def read_molecule(self, mol_filename: str):
        mol = Chem.MolFromPDBFile(mol_filename, removeHs=False)
        mol_df = mol_to_dataframe(mol)
        self.deemian_data.molecule_dataframe[mol_filename] = mol_df

    def assign_selection(self, name: str, selection: list[tuple], mol_filename: str):
        mol_df = self.deemian_data.molecule_dataframe[mol_filename]
        self.deemian_data.molecule_dataframe[name] = mol_dataframe_selection(selection, mol_df)
        selection.insert(0, mol_filename)
        self.deemian_data.selections[name] = selection

    def correct_bond(self, name: str, template: str):
        mol_df = self.deemian_data.molecule_dataframe[name]
        mol_pdb_block = dataframe_to_pdb_block(mol_df)
        self.deemian_data.molecule_pdb_block[name] = mol_pdb_block

        mol = Chem.MolFromPDBBlock(mol_pdb_block)
        template_mol = Chem.MolFromSmiles(template)
        mol = Chem.AssignBondOrdersFromTemplate(mol, template_mol)
        self.deemian_data.molecule_dataframe[name] = mol_to_dataframe(mol)

    def set_interactions(self, interactions: list[str]):
        self.deemian_data.interactions = interactions

    def set_ionizable(self, charge: str, boolean: str):
        if boolean == "true":
            boolean = True
        elif boolean == "false":
            boolean = False
        self.deemian_data.ionizable[charge] = boolean

    def set_interacting_subjects(self, subject_1: str, subject_2: str, name: str):
        self.deemian_data.interacting_subjects[name] = (subject_1, subject_2)

    def set_conformation(self, number: str):
        self.deemian_data.conformation = [int(number)]

    def set_conformation_range(self, start: str, end: str):
        self.deemian_data.conformation = list(range(int(start), int(end) + 1))

    def calculate_interactions(self, measurement_identifier: str):
        return measurement_identifier

    def write_readable_output(self, out_file: str, presentation_identifier: str):
        return (out_file, presentation_identifier)

    def write_deemian_data(self, out_file: str, presentation_identifier: str):
        return (out_file, presentation_identifier)

    def generate_deemian_data(self):
        return self.deemian_data
