from dataclasses import dataclass, field

import pandas as pd
from rdkit.Chem import AllChem as Chem
from scipy.spatial import KDTree

from deemian.chem.interaction_utility import filter_charge, generate_pair_info


@dataclass
class InteractionData:
    subject_1_mol: Chem.rdchem.Mol
    subject_2_mol: Chem.rdchem.Mol
    subject_1_df: pd.DataFrame
    subject_2_df: pd.DataFrame
    conformation_range: range
    dataframe: pd.DataFrame = field(default_factory=lambda: pd.DataFrame())

    def calculate_electrostatic(self, positive: bool, negative: bool):
        cation_mode = "all_cation" if positive else "apparent_cation"
        anion_mode = "all_anion" if negative else "apparent_anion"

        cation_1_ids = filter_charge(self.subject_1_mol, cation_mode)
        anion_1_ids = filter_charge(self.subject_1_mol, anion_mode)
        cation_2_ids = filter_charge(self.subject_2_mol, cation_mode)
        anion_2_ids = filter_charge(self.subject_2_mol, anion_mode)

        df1 = self.subject_1_df
        df2 = self.subject_2_df
        exclude_atom = ["C", "N", "S", "P"]
        cation_1_df = df1[df1.index.isin(cation_1_ids)]
        anion_1_df = df1[df1.index.isin(anion_1_ids) & ~df1["atom_symbol"].isin(exclude_atom)]
        cation_2_df = df2[df2.index.isin(cation_2_ids)]
        anion_2_df = df2[df2.index.isin(anion_2_ids) & ~df2["atom_symbol"].isin(exclude_atom)]

        for conf_num in self.conformation_range:
            conformation_column = "conf_" + str(conf_num)

            if (not cation_1_df.empty) and (not anion_2_df.empty):
                cation_1_tree = KDTree(cation_1_df[conformation_column].to_list())
                anion_2_tree = KDTree(anion_2_df[conformation_column].to_list())

                s1_as_cation = cation_1_tree.query_ball_tree(anion_2_tree, 4.5)

                pair_info_df = generate_pair_info(
                    s1_as_cation, cation_1_df, anion_2_df, conf_num, "electrostatic_cation"
                )

                if not pair_info_df.empty:
                    self.dataframe = pd.concat([self.dataframe, pair_info_df])

            if (not cation_2_df.empty) and (not anion_1_df.empty):
                cation_2_tree = KDTree(cation_2_df[conformation_column].to_list())
                anion_1_tree = KDTree(anion_1_df[conformation_column].to_list())

                s1_as_anion = anion_1_tree.query_ball_tree(cation_2_tree, 4.5)

                pair_info_df = generate_pair_info(
                    s1_as_anion, anion_1_df, cation_2_df, conf_num, "electrostatic_anion"
                )

                if not pair_info_df.empty:
                    self.dataframe = pd.concat([self.dataframe, pair_info_df])
