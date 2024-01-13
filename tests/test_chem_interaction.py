import pytest

from deemian.chem.interactions import InteractionData


@pytest.fixture
def interaction_data_n1(oseltamivir_corrected, n1, oseltamivir_corrected_df, n1_protein_A):
    interaction_data = InteractionData(oseltamivir_corrected, n1, oseltamivir_corrected_df, n1_protein_A, range(1, 2))

    return interaction_data


def test_interaction_data_calculate_electrostatic(interaction_data_n1):
    interaction_data_n1.calculate_electrostatic(positive=True, negative=True)

    df = interaction_data_n1.dataframe
    electrostatic_as_cat_df = df[df["interaction_type"] == "electrostatic_cation"]
    electrostatic_as_an_df = df[df["interaction_type"] == "electrostatic_anion"]

    assert electrostatic_as_cat_df.shape == (3, 17)
    assert electrostatic_as_an_df.shape == (12, 17)


def test_interaction_data_calculate_electrostatic_apparent_only(interaction_data_n1):
    interaction_data_n1.calculate_electrostatic(positive=False, negative=False)

    assert interaction_data_n1.dataframe.shape == (0, 0)
