from collections import namedtuple
from unittest.mock import Mock, call

from lark import Token, Tree
import pytest

from deemian.engine.director import director


Selection = namedtuple("Selection", "name selection type", defaults=["selection"])
BondCorrection = namedtuple("BondCorrection", "name template type", defaults=["bond_correction"])
ReadableOutput = namedtuple("ReadableOutput", "results out_file type", defaults=["readable_output"])
DeemianData = namedtuple("DeemianData", "out_file type", defaults=["deemian_data"])
InteractionList = namedtuple("InteractionList", "interactions type", defaults=["interaction_list"])
IncludeIonizable = namedtuple("IncludeIonizable", "charge boolean type", defaults=["include_ionizable"])
Conformation = namedtuple("Conformation", "number type", defaults=["conformation"])
ConformationRange = namedtuple("ConformationRange", "start end type", defaults=["conformation_range"])
InteractingSubject = namedtuple(
    "InteractingSubject", "subject_1 subject_2 name type", defaults=["interacting_subject"]
)


@pytest.fixture
def steps_simple() -> Tree:
    tree = Tree(
        Token("RULE", "start"),
        [
            Tree(
                "molecule_op",
                [
                    Token("IDENTIFIER", "5nzn.pdb"),
                    Selection(name="protein_A", selection=[("chain", "A"), ("and", "protein")], type="selection"),
                    Selection(
                        name="oseltamivir", selection=[("chain", "A"), ("and", "resname", "G39")], type="selection"
                    ),
                    BondCorrection(
                        name="oseltamivir",
                        template="CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)O",
                        type="bond_correction",
                    ),
                ],
            ),
            Tree(
                "measurement",
                [
                    Token("IDENTIFIER", "protein_ligand"),
                    InteractionList(interactions=["all"], type="interaction_list"),
                    IncludeIonizable(charge="positive", boolean="true", type="include_ionizable"),
                    IncludeIonizable(charge="negative", boolean="true", type="include_ionizable"),
                    InteractingSubject(
                        subject_1="protein_A",
                        subject_2="oseltamivir",
                        name="protein_A:oseltamivir",
                        type="interacting_subject",
                    ),
                    Conformation(number="1", type="conformation"),
                ],
            ),
            Tree(
                "presentation",
                [
                    Token("IDENTIFIER", "protein_ligand"),
                    ReadableOutput(results=["interactions"], out_file="protein_ligand.txt", type="readable_output"),
                    DeemianData(out_file="protein_ligand.db", type="deemian_data"),
                ],
            ),
        ],
    )

    return tree


@pytest.fixture
def steps_multiselect() -> Tree:
    tree = Tree(
        Token("RULE", "start"),
        [
            Tree(
                "molecule_op",
                [
                    Token("IDENTIFIER", "7u0n.pdb"),
                    Selection(
                        name="ace2_hotspot31",
                        selection=[("chain", "A"), ("and", "protein"), ("and", "resid_range", "1", "55")],
                        type="selection",
                    ),
                    Selection(
                        name="ace2_hotspot353",
                        selection=[("chain", "A"), ("and", "protein"), ("and", "resid_range", "341", "364")],
                        type="selection",
                    ),
                    Selection(
                        name="spike_rbm",
                        selection=[("chain", "E"), ("and", "protein"), ("and", "resid_range", "470", "510")],
                        type="selection",
                    ),
                ],
            ),
            Tree(
                "measurement",
                [
                    Token("IDENTIFIER", "ace2_spike_rbd"),
                    InteractionList(interactions=["all"], type="interaction_list"),
                    IncludeIonizable(charge="positive", boolean="true", type="include_ionizable"),
                    IncludeIonizable(charge="negative", boolean="true", type="include_ionizable"),
                    InteractingSubject(
                        subject_1="spike_rbm", subject_2="ace2_hotspot31", name="rbm_h31", type="interacting_subject"
                    ),
                    InteractingSubject(
                        subject_1="spike_rbm", subject_2="ace2_hotspot353", name="rbm_h353", type="interacting_subject"
                    ),
                    InteractingSubject(
                        subject_1="ace2_hotspot31",
                        subject_2="ace2_hotspot353",
                        name="internal_ace2",
                        type="interacting_subject",
                    ),
                    InteractingSubject(
                        subject_1="spike_rbm", subject_2="spike_rbm", name="internal_rbm", type="interacting_subject"
                    ),
                ],
            ),
            Tree(
                "presentation",
                [
                    Token("IDENTIFIER", "ace2_spike_rbd"),
                    ReadableOutput(
                        results=["interactions"], out_file="ace2_spike_rbd_detailed.txt", type="readable_output"
                    ),
                    DeemianData(out_file="ace2_spike_rbd_detailed.db", type="deemian_data"),
                ],
            ),
        ],
    )

    return tree


@pytest.fixture
def steps_multiconf() -> Tree:
    tree = Tree(
        Token("RULE", "start"),
        [
            Tree(
                "molecule_op",
                [
                    Token("IDENTIFIER", "2kw3.pdb"),
                    Selection(name="vps4", selection=[("protein",), ("and", "chain", "A")], type="selection"),
                    Selection(name="chmp6", selection=[("protein",), ("and", "chain", "B")], type="selection"),
                ],
            ),
            Tree(
                "measurement",
                [
                    Token("IDENTIFIER", "vps4_chmp6"),
                    InteractionList(interactions=["all"], type="interaction_list"),
                    IncludeIonizable(charge="positive", boolean="true", type="include_ionizable"),
                    IncludeIonizable(charge="negative", boolean="true", type="include_ionizable"),
                    InteractingSubject(
                        subject_1="vps4", subject_2="chmp6", name="vps4:chmp6", type="interacting_subject"
                    ),
                    ConformationRange(start="1", end="20", type="conformation_range"),
                ],
            ),
            Tree(
                "presentation",
                [
                    Token("IDENTIFIER", "vps4_chmp6"),
                    ReadableOutput(results=["interactions"], out_file="vps4_chmp6.txt", type="readable_output"),
                    DeemianData(out_file="vps4_chmp6.db", type="deemian_data"),
                ],
            ),
        ],
    )

    return tree


def test_engine_director_simple(steps_simple):
    data_builder = Mock()
    director(steps_simple, data_builder)

    data_builder.assert_has_calls(
        [
            call.read_molecule("5nzn.pdb"),
            call.assign_selection("protein_A", [("chain", "A"), ("and", "protein")], "5nzn.pdb"),
            call.assign_selection("oseltamivir", [("chain", "A"), ("and", "resname", "G39")], "5nzn.pdb"),
            call.correct_bond("oseltamivir", "CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)O"),
            call.set_interactions(["all"]),
            call.set_ionizable("positive", "true"),
            call.set_ionizable("negative", "true"),
            call.set_interacting_subjects("protein_A", "oseltamivir", "protein_A:oseltamivir"),
            call.set_conformation("1"),
            call.calculate_interactions("protein_ligand"),
            call.write_readable_output("protein_ligand.txt", "protein_ligand"),
            call.write_deemian_data("protein_ligand.db", "protein_ligand"),
        ]
    )


def test_engine_director_multiselect(steps_multiselect):
    data_builder = Mock()
    director(steps_multiselect, data_builder)

    data_builder.assert_has_calls(
        [
            call.read_molecule("7u0n.pdb"),
            call.assign_selection(
                "ace2_hotspot31", [("chain", "A"), ("and", "protein"), ("and", "resid_range", "1", "55")], "7u0n.pdb"
            ),
            call.assign_selection(
                "ace2_hotspot353",
                [("chain", "A"), ("and", "protein"), ("and", "resid_range", "341", "364")],
                "7u0n.pdb",
            ),
            call.assign_selection(
                "spike_rbm", [("chain", "E"), ("and", "protein"), ("and", "resid_range", "470", "510")], "7u0n.pdb"
            ),
            call.set_interactions(["all"]),
            call.set_ionizable("positive", "true"),
            call.set_ionizable("negative", "true"),
            call.set_interacting_subjects("spike_rbm", "ace2_hotspot31", "rbm_h31"),
            call.set_interacting_subjects("spike_rbm", "ace2_hotspot353", "rbm_h353"),
            call.set_interacting_subjects("ace2_hotspot31", "ace2_hotspot353", "internal_ace2"),
            call.set_interacting_subjects("spike_rbm", "spike_rbm", "internal_rbm"),
            call.calculate_interactions("ace2_spike_rbd"),
            call.write_readable_output("ace2_spike_rbd_detailed.txt", "ace2_spike_rbd"),
            call.write_deemian_data("ace2_spike_rbd_detailed.db", "ace2_spike_rbd"),
        ]
    )


def test_engine_director_multiconf(steps_multiconf):
    data_builder = Mock()
    director(steps_multiconf, data_builder)

    data_builder.assert_has_calls(
        [
            call.read_molecule("2kw3.pdb"),
            call.assign_selection("vps4", [("protein",), ("and", "chain", "A")], "2kw3.pdb"),
            call.assign_selection("chmp6", [("protein",), ("and", "chain", "B")], "2kw3.pdb"),
            call.set_interactions(["all"]),
            call.set_ionizable("positive", "true"),
            call.set_ionizable("negative", "true"),
            call.set_interacting_subjects("vps4", "chmp6", "vps4:chmp6"),
            call.set_conformation_range("1", "20"),
            call.calculate_interactions("vps4_chmp6"),
            call.write_readable_output("vps4_chmp6.txt", "vps4_chmp6"),
            call.write_deemian_data("vps4_chmp6.db", "vps4_chmp6"),
        ]
    )