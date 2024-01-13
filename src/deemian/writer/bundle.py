from copy import deepcopy
from dataclasses import asdict
from datetime import datetime
import io
import json
import tarfile
import time


def dict_factory(x):
    exclude_fields = ("start_time", "end_time")
    return {k: v for (k, v) in x if ((v is not None) and (k not in exclude_fields))}


def write_metadata(metadata, measurement_data):
    metadata_copy = deepcopy(metadata)
    metadata_copy.measurement["interacting_subjects"] = []

    for name, subjects in measurement_data.interacting_subjects.items():
        parquet_result = name.replace(":", "_") + ".parquet"
        interacting_subject = dict(name=name, subjects=subjects, results=parquet_result)
        metadata_copy.measurement["interacting_subjects"].append(interacting_subject)

    conf_range = measurement_data.conformation_range
    metadata_copy.measurement["conformation_range"] = [conf_range.start, conf_range.stop - 1]

    # setup the original metadata, because the running time is used by main script
    metadata.end_time = time.time()
    metadata.running_time = f"{metadata.end_time - metadata.start_time:.3f} seconds"
    metadata_copy.running_time = metadata.running_time
    metadata_copy.creation_time = datetime.now().strftime("%m/%d/%Y, %H:%M:%S")

    metadata_json = json.dumps(asdict(metadata_copy, dict_factory=dict_factory), indent=2).encode("utf-8")
    metadata_io = io.BytesIO(metadata_json)
    size = len(metadata_json)

    return ("deemian.json", metadata_io, size)


def write_corrected_molecule(name, pdb_block: str):
    name = name + "_corrected.pdb"
    pdb_block = pdb_block.encode("utf-8")
    pdb_io = io.BytesIO(pdb_block)
    size = len(pdb_block)

    return (name, pdb_io, size)


def write_calculation_result(name, interaction_df):
    parquet_result = name.replace(":", "_") + ".parquet"
    parquet_io = io.BytesIO()
    interaction_df.to_parquet(parquet_io)
    size = parquet_io.getbuffer().nbytes

    return (parquet_result, parquet_io, size)


def write_bundle(manifest, out_file):
    with tarfile.open(out_file, mode="w:gz") as dd_file:
        for file in manifest:
            file_name, file_object, size = file
            tarinfo = tarfile.TarInfo(file_name)
            tarinfo.size = size
            file_object.seek(0)
            dd_file.addfile(tarinfo, file_object)
        print(f"       ... wrote Deemian data: {out_file}")
