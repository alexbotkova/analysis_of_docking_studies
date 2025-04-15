from setuptools import setup, find_packages

setup(
    name="analysis_of_docking_studies",             
    version="1.0.0",                      
    author="Alexandra Botkova",
    author_email="botkova@natur.cuni.cz",
    description="A CLI tool for analyzing protein–ligand docking results.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    packages=find_packages(),             
    install_requires=[
        "click>=8.1.8",
        "freesasa>=2.2.1",
        "numpy>=2.2.4",
        "pandas>=2.2.3",
        "pytest>=8.3.5",
        "rdkit>=2024.3.5",
    ],
    python_requires='>=3.8',

    extras_require={
    "viz": [
        "pymol>=3.0.0", # PyMOL is not available via pip; install it via conda: conda install pymol-open-source
    ],
    }
)
