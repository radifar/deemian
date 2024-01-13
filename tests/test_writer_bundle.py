from dataclasses import asdict
import io
import json
import tarfile
from unittest.mock import call, patch, MagicMock

import pytest

from deemian.engine.builder import Measurement, Metadata
from deemian.writer.bundle import (
    dict_factory,
    write_metadata,
    write_corrected_molecule,
    write_calculation_result,
    write_bundle,
)


@pytest.fixture
def metadata_obj():
    metadata = Metadata()

    return metadata


@pytest.fixture
def measurement_data():
    measurement = Measurement()
    measurement.interactions.extend(["all"])
    measurement.set_ionizable("positive", "true")
    measurement.set_ionizable("negative", "true")
    measurement.interacting_subjects["oseltamivir:protein_A"] = ("oseltamivir", "protein_A")

    return measurement


def test_writer_bundle_dict_factory(metadata_obj):
    metadata_dict = asdict(metadata_obj, dict_factory=dict_factory)
    metadata_keys = list(metadata_dict.keys())

    expected_keys = "deemian_version running_time creation_time selection measurement".split()

    assert metadata_keys == expected_keys


def test_writer_bundle_write_metadata_no_conformation(metadata_obj, measurement_data):
    metadata_file = write_metadata(metadata_obj, measurement_data)
    metadata_file_name = metadata_file[0]
    metadata_io = metadata_file[1]
    metadata_json = json.loads(metadata_io.read())

    assert metadata_file_name == "deemian.json"
    assert metadata_json["measurement"]["conformation_range"] == [1, 1]


def test_writer_bundle_write_metadata_conformation_range(metadata_obj, measurement_data):
    measurement_data.conformation_range = range(1, 21)

    metadata_file = write_metadata(metadata_obj, measurement_data)
    metadata_file_name = metadata_file[0]
    metadata_io = metadata_file[1]
    metadata_json = json.loads(metadata_io.read())

    assert metadata_file_name == "deemian.json"
    assert metadata_json["measurement"]["conformation_range"] == [1, 20]


@patch.object(io, "BytesIO")
def test_writer_bundle_write_corrected_molecule(bytesio):
    name = "oseltamivir"
    pdb_block = "HETATM    1 O1A  G39 A 503"
    size = len(pdb_block)

    corrected_molecule = write_corrected_molecule(name, pdb_block)

    bytesio.assert_called_once_with(pdb_block.encode("utf-8"))
    assert corrected_molecule[0] == "oseltamivir_corrected.pdb"
    assert corrected_molecule[1] == bytesio()
    assert corrected_molecule[2] == size


@patch.object(io, "BytesIO")
def test_writer_bundle_write_calculation_result(bytesio):
    name = "oseltamivir:protein_A"
    interaction_df = MagicMock()
    bytesio.return_value.getbuffer.return_value.nbytes = 10

    calculation_result = write_calculation_result(name, interaction_df)

    assert calculation_result[0] == "oseltamivir_protein_A.parquet"
    assert calculation_result[1] == bytesio()
    assert calculation_result[2] == 10


@patch.object(tarfile, "TarInfo")
@patch("deemian.writer.bundle.tarfile.open")
def test_writer_bundle_write_bundle(mock_open, tarinfo):
    mock_addfile = MagicMock()
    mock_open.return_value.__enter__.return_value.addfile = mock_addfile

    file_name1 = "deemian.json"
    file_object1 = MagicMock()
    size1 = 50

    file_name2 = "oseltamivir_protein_A.parquet"
    file_object2 = MagicMock()
    size2 = 100

    manifest = [(file_name1, file_object1, size1), (file_name2, file_object2, size2)]
    out_file = "n1_result.dd"

    write_bundle(manifest, out_file)

    tarinfo.assert_has_calls([call("deemian.json"), call("oseltamivir_protein_A.parquet")])

    mock_addfile.assert_has_calls([call(tarinfo(), file_object1), call(tarinfo(), file_object2)])
