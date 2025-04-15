# Changelog

The project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-04-15

### Added

- Command-line tool for analyzing protein–ligand docking results.
- Support for input via PDB code or structure file.
- PyMOL integration for visualization of ligand poses and pockets.
- Display of ligand name, docking score, and pocket ID for each pose.
- Pocket analysis with FreeSASA for SASA, GRAVY index, and net charge.
- Parsing of additional `predictions` and `residues` input files.
- `--res` flag to display only residues involved in hydrogen bonding per pose.
- Modular design enabling future extension of commands and options.

### Removed

### Fixed