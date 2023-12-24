import pytest

from deemian.chem.interactions import InteractionData


@pytest.fixture
def interaction_data_n1(oseltamivir_corrected, n1, oseltamivir_corrected_df, n1_protein_A):
    interaction_data = InteractionData(oseltamivir_corrected, n1, oseltamivir_corrected_df, n1_protein_A, [1])

    return interaction_data


def test_interaction_data_calculate_electrostatic(interaction_data_n1):
    interaction_data_n1.calculate_electrostatic(positive=True, negative=True)

    assert interaction_data_n1.electrostatic_s1_as_cation.shape == (3, 17)
    interaction_data_n1.electrostatic_s1_as_anion.shape == (12, 17)


def test_interaction_data_calculate_electrostatic_apparent_only(interaction_data_n1):
    interaction_data_n1.calculate_electrostatic(positive=False, negative=False)

    assert interaction_data_n1.electrostatic_s1_as_cation.shape == (0, 0)
    assert interaction_data_n1.electrostatic_s1_as_anion.shape == (0, 0)
