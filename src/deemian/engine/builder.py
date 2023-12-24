from collections import defaultdict
from dataclasses import dataclass, field

import pandas as pd
from rdkit.Chem import AllChem as Chem

from deemian.chem.interactions import InteractionData
from deemian.chem.reader import mol_to_dataframe
from deemian.chem.selection import mol_dataframe_selection
from deemian.chem.utility import dataframe_to_pdb_block


@dataclass
class Molecule:
    rdkit_mol: Chem.rdchem.Mol
    mol_dataframe: pd.DataFrame = None


@dataclass
class Selection:
    mol_parent: str
    mol_dataframe: pd.DataFrame
    mol_pdb_block: str = None


@dataclass
class Measurement:
    interactions: list = field(default_factory=lambda: [])
    ionizable: dict = field(default_factory=lambda: {"positive": False, "negative": False})
    interacting_subjects: dict = field(default_factory=lambda: {})
    conformation: list = field(default_factory=lambda: [])

    def conformation_range(self, start, end):
        self.conformation = list(range(int(start), int(end) + 1))

    def set_ionizable(self, charge: str, boolean: str):
        if boolean == "true":
            boolean = True
        elif boolean == "false":
            boolean = False
        self.ionizable[charge] = boolean


@dataclass
class DeemianData:
    molecule: dict[str, Molecule] = field(default_factory=lambda: {})
    selection: dict[str, Selection] = field(default_factory=lambda: {})
    measurement: dict[str, Measurement] = field(default_factory=lambda: defaultdict(Measurement))
    interaction_details: dict = field(default_factory=lambda: {})
    readable_output: dict = field(default_factory=lambda: {})

    def add_molecule(self, name):
        mol = Chem.MolFromPDBFile(name, removeHs=False)
        mol_df = mol_to_dataframe(mol)
        self.molecule[name] = Molecule(mol, mol_df)

    def add_selection(self, name, selection, mol_parent):
        parent_df = self.molecule[mol_parent].mol_dataframe
        selection_df = mol_dataframe_selection(selection, parent_df)

        self.selection[name] = Selection(mol_parent, selection_df)

    def correct_bond(self, name, template):
        selection_df = self.selection[name].mol_dataframe
        selection_pdb_block = dataframe_to_pdb_block(selection_df)
        selection_mol = Chem.MolFromPDBBlock(selection_pdb_block)

        template_mol = Chem.MolFromSmiles(template)
        corrected_mol = Chem.AssignBondOrdersFromTemplate(template_mol, selection_mol)
        corrected_df = mol_to_dataframe(corrected_mol)

        self.molecule[name] = Molecule(corrected_mol)
        self.selection[name] = Selection(name, corrected_df, selection_pdb_block)

    def add_measurement(self, id):
        return self.measurement[id]

    def calculate_interactions(self, id):
        measurement = self.measurement[id]

        for pair in measurement.interacting_subjects:
            subject_1, subject_2 = measurement.interacting_subjects[pair]
            subject_1 = self.selection[subject_1]
            subject_2 = self.selection[subject_2]
            subject_1_mol = self.molecule[subject_1.mol_parent]
            subject_2_mol = self.molecule[subject_2.mol_parent]
            subject_1_df = subject_1.mol_dataframe
            subject_2_df = subject_2.mol_dataframe
            conformation = measurement.conformation
            if not any(conformation):
                conformation = [1]

            interaction_data = InteractionData(subject_1_mol, subject_2_mol, subject_1_df, subject_2_df, conformation)

            interaction_type = measurement.interactions

            if ("electrostatic" in interaction_type) or ("all" in interaction_type):
                interaction_data.calculate_electrostatic(**measurement.ionizable)

            self.interaction_details[pair] = interaction_data

    def write_readable_output(self, out_file: str, presentation_id: str):
        return (out_file, presentation_id)

    def write_deemian_data(self, out_file: str, presentation_id: str):
        return (out_file, presentation_id)
