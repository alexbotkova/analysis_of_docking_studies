import pytest
import re
import pandas as pd
from io import StringIO
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../dockinspect')))
from pymol_tools import get_pocket_residues_dict, get_pocket_atomids_dict
from my_parser import get_df

predictions_filepath = "/Users/alexbotkova/analysis_of_docking_studies/test_files/urea/prankweb-2SRC/structure.cif_predictions.csv"

@pytest.fixture
def predictions_df():
    return get_df(predictions_filepath)

MOCK_PREDICTIONS_WITH_ATOMS = StringIO("""name,residue_ids,surf_atom_ids
pocket1,A_123 A_124 A_125,101 102 103
pocket2,B_200 B_201,201 202
""")

def test_get_pocket_residues_dict(predictions_df):
    pocket_residues = get_pocket_residues_dict(predictions_df)

    assert isinstance(pocket_residues, dict)
    assert len(pocket_residues) > 0

    for pocket, selection in pocket_residues.items():
        assert isinstance(pocket, str)
        assert isinstance(selection, str)
        assert "chain" in selection and "resi" in selection
        assert "+" in selection or re.search(r'resi \d+', selection)

    expected_plus_counts = {
        'pocket1': 37,
        'pocket2': 13,
        'pocket3': 7,
        'pocket4': 9,
        'pocket5': 6,
        'pocket6': 8
    }

    actual_plus_counts = {
        pocket: selection.count("+")
        for pocket, selection in pocket_residues.items()
        if pocket in expected_plus_counts
    }

    assert actual_plus_counts == expected_plus_counts, f"Mismatch in '+' counts: {actual_plus_counts}"

def test_get_pocket_atomids_dict(predictions_df):
    pocket_atoms = get_pocket_atomids_dict(predictions_df)

    assert isinstance(pocket_atoms, dict)
    assert len(pocket_atoms) > 0

    for pocket, selection in pocket_atoms.items():
        assert isinstance(pocket, str)
        assert isinstance(selection, str)
        assert selection.startswith("id ")
        assert all(x.isdigit() for x in selection.replace("id ", "").split("+"))

def test_get_pocket_residues_dict_manual():
    MOCK_PREDICTIONS_WITH_ATOMS.seek(0)
    df = pd.read_csv(MOCK_PREDICTIONS_WITH_ATOMS)

    result = get_pocket_residues_dict(df)
    print(result)
    expected = {
        "pocket1": "chain A and resi 123+124+125",
        "pocket2": "chain B and resi 200+201"
    }
    assert result == expected

def test_get_pocket_atomids_dict_manual():
    MOCK_PREDICTIONS_WITH_ATOMS.seek(0)
    df = pd.read_csv(MOCK_PREDICTIONS_WITH_ATOMS)

    result = get_pocket_atomids_dict(df)
    print(result)
    expected = {
        "pocket1": "id 101+102+103",
        "pocket2": "id 201+202"
    }
    assert result == expected