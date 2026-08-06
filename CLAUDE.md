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
  → `(keep: List[int], delete: List[int], replace: List[int])`. Builds a
  connectivity matrix from bonds (guessing if absent), classifies
  heavy/sp3/oxygen atoms, walks chains from terminal atoms, returns atom
  index lists.
- **`crop_chains(universe, chain_type='alkyl', max_chain_length=1)`
  returns a 3-tuple `(new_universe, keep, replace)`, not a bare
  `mda.Universe`** — its docstring originally only documented the
  `mda.Universe` return, which is exactly what caused the
  `crystal_analyzer` tuple-unpacking bug this session (its `__init__`
  and `export_dimers()` both originally assumed a bare `Universe`; see
  `crystal_analyzer`'s own CLAUDE.md). **Fixed 2026-08-06:** docstring
  and return type hint corrected to document the actual 3-tuple return.

`TrajectoryProcessor(ChainCropper)` — reuses cropping indices computed
once from the topology across all trajectory frames;
`process_trajectory(structure_file, output_path, trajectory_file=None,
chain_type=, max_chain_length=, n_jobs=-1)` is the main entry point used
by both the simple Python API and the CLI. Serial or joblib-parallel
per-frame processing with tqdm progress bars.

## Input format

Structure (via MDAnalysis): GRO, PDB, MOL2, XYZ. Trajectory: XTC, TRR,
DCD.

## Output format

Structure: GRO, PDB, XYZ. Trajectory: XTC, TRR, DCD — chosen by the
output path's extension.

## Known gaps / TODOs

- **Likely bug (not executed/confirmed this session), serial trajectory
  path:** `TrajectoryProcessor._apply_cropping()` (line ~473) calls
  `universe.add_TopologyAttr('name', names)` /
  `add_TopologyAttr('element', elements)` unconditionally every time it
  runs. `_write_trajectory()`'s serial branch (`n_jobs == 1`, line
  ~626-632) calls `_apply_cropping(universe)` once per frame on the
  *same* `universe` object — MDAnalysis raises `ValueError` when adding a
  `TopologyAttr` that already exists, so this would plausibly crash on
  the second frame of any multi-frame trajectory processed serially.
  The parallel branch (`n_jobs != 1`) avoids this because it only calls
  `_apply_cropping` once (for the first frame) and processes the rest via
  the standalone `_process_frame()`, which doesn't touch `TopologyAttr`s.
  Worth a real end-to-end test with a multi-frame trajectory + `-j 1`
  before trusting the serial path.
- `chain_cropper`'s own `identify_chains_to_crop()` is actually the
  correct, working implementation of `ether`-aware chain identification
  (it genuinely includes oxygen atoms in the walk, unlike the similarly-named
  `get_sp2(..., ether=True)` helpers duplicated in `SelIntCoords` and
  `oligomer_builder`, which silently no-op on that branch) — if those
  packages' ether support is ever needed for real, this module is the
  reference implementation to port from, not the other way round.
- **Fixed 2026-08-06:** `setup.py`'s `install_requires` listed only
  `MDAnalysis`/`numpy`, but `tqdm` and `joblib` are hard runtime imports
  (`_write_trajectory`, `process_trajectory`) — added both.

## Dependencies

MDAnalysis>=2.0.0, numpy>=1.19.0 (declared); `tqdm`, `joblib` (used at
runtime, undeclared — see above). No vendored code in this repo.
