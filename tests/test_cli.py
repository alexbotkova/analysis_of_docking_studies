import pytest
from click.testing import CliRunner
from dockinspect.main import cli

runner = CliRunner()

@pytest.mark.parametrize("args, expected_keywords", [
    (["ligand_info", "NC(=O)N"], ["logP", "SASA", "TPSA", "Volume", "Charge"]),
    (["pocket_info", "--help"], ["SASA", "GRAVY", "charge"]),
    (["poses_info", "--help"], ["--csv", "--res"]),
    (["visualize", "--help"], ["PyMOL", "pose", "mode", "distance", "angle"]),
    (["shell", "--help"], ["interactive shell", "ligand_smiles", "pdb_code", "vina_file"]),
])
def test_cli_commands_show_keywords(args, expected_keywords):
    result = runner.invoke(cli, args)
    assert result.exit_code == 0
    for keyword in expected_keywords:
        assert keyword in result.output

def test_pocket_info_missing_files():
    """Test pocket_info fails gracefully with missing files."""
    result = runner.invoke(cli, ["pocket_info", "missing.pdbqt", "missing.csv", "missing.csv"])
    assert result.exit_code != 0
    assert "Error" in result.output or "No such file" in result.output

def test_shell_command_no_input():
    """Test that shell with missing arguments fails properly."""
    result = runner.invoke(cli, ["shell"])
    assert result.exit_code != 0
    assert "Missing argument" in result.output

def test_visualize_command_missing_files():
    """Test visualize gracefully handles missing inputs."""
    result = runner.invoke(cli, ["visualize", "missing", "missing.pdbqt"])
    assert result.exit_code != 0
    assert "Error" in result.output or "No such file" in result.output
