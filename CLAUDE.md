@../pipeline-docs/CLAUDE.md

# CLAUDE.md

## What this repo is

Removes/truncates alkyl or ether side chains from a molecular structure
(single frame or full MD trajectory), capping the resulting broken bonds
with hydrogen at a configurable distance. Prepares structures for
electronic-structure/coupling calculations where long aliphatic side
chains are irrelevant.

## Role in the pipeline

Step 3 — trajectory postprocessing (general case; `crystal_analyzer`
handles the crystal-specific case and imports this package directly).
Upstream: an equilibrated MD trajectory/structure from step 2. Downstream:
`pyscf_tints` (step 4, via a dimer/monomer-pair extraction step —
`oligomer_builder`'s breaker or manual selection) or `crystal_analyzer`.

## Public interface

Console script `chain-cropper` (`chain_cropper.cli:main`):
`chain-cropper topology [trajectory] [-o OUT] [--chain-type alkyl|ether]
[--max-length N] [--cap-distance 1.09] [-j N_JOBS] [--batch] [-v]`.
Structure-only or structure+trajectory modes; `--batch` processes all
topology files found in the current directory. Output auto-named
`<stem>_cropped<ext>` unless `-o` given.

Python API — `ChainCropper(cap_distance=1.09)`:
- `identify_chains_to_crop(universe, chain_type='alkyl', max_chain_length=1)`
  → `(keep, delete, replace)` atom-index lists. Builds a connectivity
  matrix from bonds (guessing if absent), classifies heavy/sp3/oxygen
  atoms, walks chains from terminal atoms.
- **`crop_chains(...)` returns `(new_universe, keep, replace)`, not a
  bare `mda.Universe`** — check the return arity before unpacking.

`TrajectoryProcessor(ChainCropper)` — reuses cropping indices computed
once from the topology across all trajectory frames.
`process_trajectory(structure_file, output_path, trajectory_file=None,
chain_type=, max_chain_length=, n_jobs=-1, structure_universe=None)` is
the main entry point, used by both the Python API and the CLI.
- `structure_universe`: pass a Universe carrying real bonds (e.g. from a
  GROMACS topology) so side-chain identification uses true connectivity
  instead of guessing it off geometry — both more reliable and much
  faster on a large system. `structure_file` must still describe the
  same atoms in the same order (used for coordinates/pairing).
- `final_indices`: after `_apply_cropping` runs, the ORIGINAL-structure
  index of each atom that survives into the cropped structure, in order
  — so cropped atom `i` came from original atom `final_indices[i]`.
  Exact because cropping keeps atoms in sorted original-index order and
  *repurposes* a deleted heavy atom as each capping hydrogen rather than
  appending a new one. Lets a caller (`polymer_couplings`) remap
  per-atom data (e.g. real bonds) across the crop without reading it
  back off geometry.

## Input/output format

Structure (MDAnalysis): GRO, PDB, MOL2, XYZ. Trajectory: XTC, TRR, DCD —
chosen by the output path's extension.

## Known gaps

- `identify_chains_to_crop`'s own `ether`-aware chain walk is the
  correct, working implementation — the similarly-named `get_sp2(...,
  ether=True)` helpers duplicated in `SelIntCoords`/`oligomer_builder`
  silently no-op on that branch. Port from here if their ether support
  is ever needed for real.
- No CI; correctness has been verified by live-testing against real MD
  trajectories rather than a full test suite (see `tests/` for what
  coverage does exist — currently the cropping-logic regression test).

## Dependencies

MDAnalysis>=2.0.0, numpy>=1.19.0, tqdm, joblib. No vendored code.
