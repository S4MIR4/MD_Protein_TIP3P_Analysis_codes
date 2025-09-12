# MD_Protein_TIP3P_Analysis_codes
# Protein MD Simulation Analysis in TIP3P Water

This repository provides Python codes for analyzing molecular dynamics (MD) simulations of proteins solvated in TIP3P water. The scripts extract structural, dynamic, and energetic insights into protein–solvent interactions, enabling comprehensive biophysical characterization.

## Features
- **Structural Analyses**
  - RMSD & RMSF: Global stability and residue flexibility
  - Radius of Gyration (Rg): Protein compactness
  - Secondary Structure Evolution (via DSSP)
  - Distance, angle, and dihedral analyses
- **Protein–Water Interactions**
  - Hydrogen bond counts and lifetimes
  - Radial distribution functions (RDFs)
  - Solvent accessible surface area (SASA)
  - Water residence times and hydration shell structure
- **Energetics & Electrostatics**
  - Dipole moment calculations and time evolution
  - Free energy landscapes (FEL) along RMSD, Rg, or PCA
  - Dielectric constant estimation from dipole fluctuations
- **Correlation & Dynamics**
  - Principal Component Analysis (PCA) of collective motions
  - Dynamic Cross-Correlation Maps (DCCM)
  - Bond distance and angular correlations with symmetry
- **Visualization**
  - Publication-quality plots for all major analyses
  - 2D/3D free energy surfaces and density maps

This repository will also provide the backbone structure for MD simulation on Gromacs and calculation of aforementioned properties
molecular-dynamics-repo/
│
├── README.md                # Overview of your project
├── setup/                   # Installation, input preparation
│   ├── system_setup.md
│   ├── topology_generation.md
│   └── solvation_ions.md
│
├── simulation/              # Command lines for running MD
│   ├── energy_minimization.md
│   ├── nvt_equilibration.md
│   ├── npt_equilibration.md
│   └── production_run.md
│
├── analysis/                # Post-simulation calculations
│   ├── rmsd.md
│   ├── rmsf.md
│   ├── gyration.md
│   ├── hydrogen_bonds.md
│   └── energy_profile.md
│
└── utilities/               # General commands / helper scripts
    ├── file_conversion.md
    ├── trajectory_processing.md
    └── visualization.md
