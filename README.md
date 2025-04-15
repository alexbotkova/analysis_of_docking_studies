# DockInspect

**DockInspect** is a Python tool for analyzing and visualizing protein–ligand docking results.  
It provides command-line utilities and visualization features built on top of PyMOL and FreeSASA.

---

## Features

- Extract and summarize docking poses
- Compute pocket properties like SASA, GRAVY, and net charge
- Visualize docking results with PyMOL
- Command-line interface for quick access
- Modular design for easy integration

---

## Installation

> Requires Python 3.8+

```bash
git clone https://github.com/yourusername/dockinspect.git
cd dockinspect
# For editable development mode: pip install -e .
pip install . 
``` 

## Tests

Tests can be run using
```bash
pytest tests
``` 