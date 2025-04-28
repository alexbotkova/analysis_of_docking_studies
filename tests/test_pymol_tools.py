import pytest
import pandas as pd
from io import StringIO
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../dockinspect')))
from pymol_tools import get_pocket_residues_dict
from my_csv_parser import get_df

MOCK_PREDICTIONS_WITH_ATOMS = StringIO("""name,residue_ids,surf_atom_ids
pocket1,A_123 A_124 A_125,101 102 103
pocket2,B_200 B_201,201 202
""")

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