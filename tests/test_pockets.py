import pytest
from unittest.mock import patch

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../dockinspect')))
from pockets import Pockets

@pytest.fixture
def mock_pocket_locations():
    return {
        "pocket1": ["A_123", "A_124", "A_125"],
        "pocket2": ["B_200", "B_201"]
    }

@pytest.fixture
def mock_location_aa():
    return {
        "A_123": "ARG",  # +1
        "A_124": "GLY",  #  0
        "A_125": "ASP",  # -1
        "B_200": "HIS",  # +1
        "B_201": "LYS"   # +1
    }

def test_get_pocket_charge(mock_pocket_locations, mock_location_aa):
    charges = Pockets.get_pocket_charge(mock_pocket_locations, mock_location_aa)
    assert charges == {
        "pocket1": 0,  # +1 +0 -1
        "pocket2": 2   # +1 +1
    }

def test_get_pocket_gravy(mock_pocket_locations, mock_location_aa):
    gravys = Pockets.get_pocket_gravy(mock_pocket_locations, mock_location_aa)
    expected_gravys = {
        "pocket1": round((-4.5 + -0.4 + -3.5) / 3, 2),  # ARG, GLY, ASP
        "pocket2": round((-3.2 + -3.9) / 2, 2)          # HIS, LYS
    }
    assert gravys == expected_gravys