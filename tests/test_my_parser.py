import os
import pandas as pd
import pytest
from io import StringIO

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../dockinspect')))
from my_parser import get_df, parse_predictions, parse_residues

predictions_filepath = "/Users/alexbotkova/analysis_of_docking_studies/test_files/urea/prankweb-2SRC/structure.cif_predictions.csv"
residues_filepath = "/Users/alexbotkova/analysis_of_docking_studies/test_files/urea/prankweb-2SRC/structure.cif_residues.csv"

MOCK_PREDICTIONS = StringIO("""name,residue_ids
pocket1,A_123 A_124 A_125
pocket2,B_200 B_201
""")

MOCK_RESIDUES = StringIO("""chain,residue_label,residue_name
A,123,ARG
A,124,GLY
A,125,ASP
B,200,HIS
B,201,LYS
""")

def test_get_df():
    df = get_df(predictions_filepath)
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert "name" in df.columns
    assert "residue_ids" in df.columns

def test_parse_predictions():
    pocket_dict = parse_predictions(predictions_filepath)
    assert isinstance(pocket_dict, dict)
    assert len(pocket_dict) > 0
    for key, val in pocket_dict.items():
        assert isinstance(key, str)
        assert isinstance(val, list)
        assert all(isinstance(r, str) for r in val)

def test_parse_residues():
    residue_dict = parse_residues(residues_filepath)
    assert isinstance(residue_dict, dict)
    assert len(residue_dict) > 0
    for key, val in residue_dict.items():
        assert "_" in key  
        assert isinstance(val, str)

def test_parse_predictions_manual(tmp_path):
    mock_path = tmp_path / "mock_preds.csv"
    df = pd.read_csv(MOCK_PREDICTIONS)
    df.columns = df.columns.str.strip()
    df.to_csv(mock_path, index=False)

    expected = {
        "pocket1": ["A_123", "A_124", "A_125"],
        "pocket2": ["B_200", "B_201"]
    }
    result = parse_predictions(mock_path)
    assert result == expected

def test_parse_residues_manual(tmp_path):
    mock_path = tmp_path / "mock_residues.csv"
    df = pd.read_csv(MOCK_RESIDUES)
    df.columns = df.columns.str.strip()
    df.to_csv(mock_path, index=False)

    expected = {
        "A_123": "ARG",
        "A_124": "GLY",
        "A_125": "ASP",
        "B_200": "HIS",
        "B_201": "LYS"
    }
    result = parse_residues(mock_path)
    assert result == expected