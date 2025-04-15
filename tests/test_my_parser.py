import os
import pandas as pd
import pytest
from io import StringIO

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../dockinspect')))
from my_parser import get_df, parse_predictions, parse_residues

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