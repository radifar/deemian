from pathlib import Path

import pytest

from deemian.engine import processor


@pytest.fixture
def parse_7qxy_i1g():
    """Parse the script for simple analysis of 7QXY"""

    base_path = Path(__file__).parent

    with open(base_path / "data/01-7qxy-i1g.txt") as reader:
        text = reader.read()
        steps = processor.parser(text)

        return steps


def test_parser_open_pdb_complex(parse_7qxy_i1g):
    """
    open 7qxy.pdb as target
    =======================

    open_molecule
    ├── 7qxy.pdb
    └── target

    Tree('open_molecule',
    [
        Token('FILE_NAME', '7qxy.pdb'),
        Token('IDENTIFIER', 'target')
    ])
    """

    step_open = parse_7qxy_i1g.children[0]

    file_name = step_open.children[0]
    target_molecule = step_open.children[1]

    assert step_open.data == "open_molecule"
    assert file_name == "7qxy.pdb"
    assert target_molecule == "target"


def test_parser_select_protein(parse_7qxy_i1g):
    """
    select prot [
        selection protein
        molecule target
    ]
    =====================

    selection
    ├── prot
    ├── selection_string
    │   └── protein
    └── target_molecule
        └── target

    Tree('selection',
    [
        Token('IDENTIFIER', 'prot'),
        Tree('selection_string',
        [
            Token('MACRO', 'protein')
        ]),
        Tree('target_molecule',
        [
            Token('IDENTIFIER', 'target')
        ])
    ])
    """

    step_select = parse_7qxy_i1g.children[1]

    selection_identifier = step_select.children[0]
    selection_string = step_select.children[1].children
    target_molecule = step_select.children[2].children

    assert step_select.data == "selection"
    assert selection_identifier == "prot"
    assert selection_string == ["protein"]
    assert target_molecule == ["target"]


def test_parser_select_ligand(parse_7qxy_i1g):
    """
    select dcp-pd3 [
        selection resname I1G
        molecule target
    ]
    =========================

    selection
    ├── dcp-pd3
    ├── selection_string
    │   └── resname
    │       └── I1G
    └── target_molecule
        └── target

    Tree('selection',
    [
        Token('IDENTIFIER', 'dcp-pd3'),
        Tree('selection_string',
        [
            Tree('resname',
            [
                Token('IDENTIFIER', 'I1G')
            ])
        ]),
        Tree('target_molecule',
        [
            Token('IDENTIFIER', 'target')
        ])
    ])
    """

    step_select = parse_7qxy_i1g.children[2]

    selection_identifier = step_select.children[0]
    selection_string = step_select.children[1].children[0]
    target_molecule = step_select.children[2].children[0]

    assert step_select.data == "selection"
    assert selection_identifier == "dcp-pd3"
    assert selection_string.data == "resname"
    assert selection_string.children == ["I1G"]
    assert target_molecule == "target"


def test_parser_measurement(parse_7qxy_i1g):
    """
    measure protein_ligand [
        interactions all
        between prot and dcp-pd3
    ]
    ============================

    measurement
    ├── protein_ligand
    ├── interaction_list
    │   └── all
    └── interacting_subject
        ├── between
        ├── prot
        ├── and
        └── dcp-pd3

    Tree('measurement',
    [
        Token('IDENTIFIER', 'protein_ligand'),
        Tree('interaction_list',
        [
            Token('INTERACTION', 'all')
        ]),
        Tree('interacting_subject',
        [
            Token('PREPOSITION', 'between'),
            Token('IDENTIFIER', 'prot'),
            Token('PREPOSITION', 'and'),
            Token('IDENTIFIER', 'dcp-pd3')
        ])
    ])
    """

    step_measurement = parse_7qxy_i1g.children[3]

    measurement_identifier = step_measurement.children[0]
    interaction_list = step_measurement.children[1].children
    interacting_subject = step_measurement.children[2].children

    assert step_measurement.data == "measurement"
    assert measurement_identifier == "protein_ligand"
    assert interaction_list == ["all"]
    assert interacting_subject == ["between", "prot", "and", "dcp-pd3"]


def test_parser_presentation(parse_7qxy_i1g):
    """
    present protein_ligand [
        interactions protein_ligand.txt
    ]
    ===================================

    presentation
    ├── protein_ligand
    └── presentation_format
        ├── interactions
        └── protein_ligand.txt

    Tree('presentation',
    [
        Token('IDENTIFIER', 'protein_ligand'),
        Tree('presentation_format',
        [
            Token('RESULTS', 'interactions'),
            Token('IDENTIFIER', 'protein_ligand.txt')
        ])
    ])
    """

    step_presentation = parse_7qxy_i1g.children[4]

    presentation_identifier = step_presentation.children[0]
    presentation_format = step_presentation.children[1]

    assert step_presentation.data == "presentation"
    assert presentation_identifier == "protein_ligand"
    assert presentation_format.data == "presentation_format"
    assert presentation_format.children == ["interactions", "protein_ligand.txt"]
