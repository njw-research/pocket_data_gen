# Enzyme Engineering Pipeline

A computational pipeline for enzyme engineering focusing on lipase design through pocket prediction, molecular docking, and molecular dynamics.

## Overview

This project implements a three-stage computational workflow for enzyme engineering:
1. Pocket identification using PocketGen
2. Binding affinity assessment with AutoDock Vina
3. Structural refinement via OpenMM molecular dynamics

## Pipeline Stages

### 1. Pocket Generation
- Input: 50 candidate enzymes structurally similar to references (6or3, 6osz, 6or8)
- Target ligands: triacylglycerol, diacylglycerol, carboxylate
- Uses PocketGen for binding site prediction
- Generates dataset of potential reactive pockets

### 2. Binding Analysis
- Molecular docking using AutoDock Vina
- Configurable protein flexibility
- Ranks pockets based on binding affinity
- Identifies promising binding conformations

### 3. MD Refinement
- OpenMM molecular dynamics simulations
- Rigid body constraints for transition state ligand
- Local conformational sampling
- Structure relaxation around binding sites
- Final ranking based on transition state stabilization

## Installation

```bash
# Create and activate conda environment
conda create -n enzyme_design python=3.9
conda activate enzyme_design

# Install dependencies using uv for speed
uv pip install -r requirements.txt
```

## Usage

[Instructions to be added as pipeline components are implemented]

## Dependencies

- PocketGen
- AutoDock Vina
- OpenMM
- Additional requirements in requirements.txt

## Project Structure

```
enzyme_pipeline/
├── src/
│   ├── pocket_prediction/
│   ├── docking/
│   └── md/
├── tests/
├── data/
│   ├── structures/
│   ├── pockets/
│   └── results/
└── configs/
```

## License

[License details to be added]

## Contributors

[Add team members]

## Citation

If you use this pipeline in your research, please cite:
[Citation details to be added]
