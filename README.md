# ChainCropper

A Python package for cropping molecular side chains from structures and trajectories using MDAnalysis.

## Overview

ChainCropper provides functionality to identify and remove alkyl or ether side chains from molecular structures while properly capping broken bonds with hydrogen atoms. It supports both single structures and molecular dynamics trajectories.

## Features

- **Configurable Chain Length**: Control how many heavy atoms to keep in side chains
- **Trajectory Processing**: Process entire MD trajectories frame by frame
- **Multiple File Formats**: Support for GRO, PDB, XYZ, XTC, TRR, DCD formats
- **Automatic Bond Capping**: Properly cap broken bonds with hydrogen atoms
- **Batch Processing**: Process multiple files at once

## Installation

### From Source

```bash
git clone <repository-url>
cd chain_cropper
pip install -e .
```

### Requirements

- Python >= 3.7
- MDAnalysis >= 2.0.0
- NumPy

## Usage

### Command Line Interface

#### Basic Usage

```bash
# Crop alkyl chains, keeping 1 carbon atom
chain-cropper input.gro -o output.gro

# Crop ether chains, keeping 2 heavy atoms
chain-cropper input.gro -o output.gro --chain-type ether --max-length 2

# Remove all side chains
chain-cropper input.pdb -o output.xyz --max-length 0
```

#### Trajectory Processing

```bash
# Process entire trajectory
chain-cropper topology.gro trajectory.xtc -o output.xtc --chain-type alkyl --max-length 1
```

#### Batch Processing

```bash
# Process all structure files in current directory
chain-cropper --batch -o output_directory/
```

#### Command Line Options

- `topology`: Input topology file (GRO, PDB, MOL2)
- `trajectory`: Optional trajectory file (XTC, TRR, DCD)
- `-o, --output`: Output file path or directory (for batch mode)
- `--chain-type`: Type of chains to crop (`alkyl` or `ether`, default: `alkyl`)
- `--max-length`: Maximum heavy atoms to keep in chains (default: 1)
- `--cap-distance`: Distance for capping hydrogens in Å (default: 1.09)
- `--batch`: Process all files in current directory
- `--verbose`: Enable verbose output

### Python API

```python
from chain_cropper import ChainCropper, TrajectoryProcessor
import MDAnalysis as mda

# Load your structure
universe = mda.Universe('input.gro')

# Initialize cropper
cropper = ChainCropper(cap_distance=1.09)

# Crop chains
cropped_universe = cropper.crop_chains(
    universe, 
    chain_type='alkyl', 
    max_chain_length=1
)

# Save result
cropped_universe.atoms.write('output.gro')

# For trajectories
processor = TrajectoryProcessor(cropper)
processor.process_trajectory(
    universe, 
    'output.xtc', 
    chain_type='alkyl', 
    max_chain_length=1
)
```

## Algorithm Details

### Chain Identification

1. **Connectivity Analysis**: Builds connectivity matrix from molecular bonds
2. **Atom Classification**: Identifies heavy atoms, sp3 carbons, and oxygen atoms
3. **Terminal Detection**: Finds chain atoms connected to non-chain (core) atoms
4. **Chain Tracing**: Follows connectivity to build complete chain sequences

### Chain Cropping Process

1. **Chain Mapping**: Identifies which atoms to keep, delete, or replace
2. **Bond Breaking**: Removes atoms beyond the specified chain length
3. **Hydrogen Capping**: Adds hydrogen atoms to cap broken bonds
4. **Structure Rebuilding**: Creates new universe with modified connectivity

## Supported File Formats

### Input Formats
- **Structure**: GRO, PDB, MOL2
- **Trajectory**: XTC, TRR, DCD

### Output Formats
- **Structure**: GRO, PDB, XYZ
- **Trajectory**: XTC, TRR, DCD

## Examples

### Example 1: Basic Alkyl Chain Cropping

```bash
# Input: Molecule with long alkyl chains
# Output: Same molecule with chains truncated to 1 carbon
chain-cropper molecule.gro -o cropped.gro --max-length 1
```

### Example 2: Ether Chain Processing

```bash
# Process ether-containing molecule, keeping 2 heavy atoms in chains
chain-cropper input.pdb -o output.pdb --chain-type ether --max-length 2
```

### Example 3: Trajectory Analysis

```bash
# Process MD trajectory, cropping chains in each frame
chain-cropper topology.gro trajectory.xtc -o cropped_traj.xtc
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Changelog

### Version 0.1
- Initial release
- Support for alkyl and ether chain cropping
- Trajectory processing capabilities
- Command line interface
- Batch processing mode
