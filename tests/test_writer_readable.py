from unittest.mock import Mock, patch, mock_open

import pandas as pd
import pytest

from deemian.chem.interactions import InteractionData
from deemian.engine.builder import Measurement
from deemian.writer.readable import generate_report, write_readable


@pytest.fixture
def results_g39_protein_A():
    return pd.read_parquet("tests/data/results_g39_protein_A.parquet.gzip")


def test_writer_readable_report(results_g39_protein_A):
    subject_1_mol = Mock()
    subject_2_mol = Mock()
    subject_1_df = pd.DataFrame()
    subject_2_df = pd.DataFrame()
    conformation = []

    g39_protein_A_interaction_data = InteractionData(
        subject_1_mol, subject_2_mol, subject_1_df, subject_2_df, conformation, results_g39_protein_A
    )

    subjects = {"g39:protein_A": ("g39", "protein_A")}
    results = {"g39:protein_A": g39_protein_A_interaction_data}
    measurement = Measurement(interacting_subjects=subjects, calculation_results=results)
    report = generate_report(measurement, "detailed_conf_first")

    assert "G39" in report
    assert "protein_A" in report


def test_writer_readable_write_readable():
    m = mock_open()
    with patch("builtins.open", m):
        content = "writing reports"
        out_file = "protein_ligand.txt"
        write_readable(content, out_file)

    handle = m()
    handle.write.assert_called_once_with("writing reports")
