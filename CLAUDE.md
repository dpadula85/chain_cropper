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
`<stem>_cropped<ext>` unless `-o` given — `-o` may be a directory (both
outputs named `<stem>_cropped<ext>` inside it) or a single filename (see
Known gaps for the file case's exact derivation rule).

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

`TrajectoryProcessor.process_trajectory(..., structure_universe=None)`
(**added 2026-08-06**) — pass a Universe carrying real bonds (from a
GROMACS topology, say) and the crop decides what counts as side chain
from the true connectivity instead of guessing it off the geometry.
Guessing both misses real bonds and invents absent ones on a strained MD
snapshot, and it is the slow step on a large system. `structure_file` is
still used for coordinates and for pairing with the trajectory, so the
two must describe the same atoms in the same order. Purely additive: the
default path is unchanged.

`TrajectoryProcessor.final_indices` (**added 2026-08-06**) — after
`_apply_cropping` has run, the ORIGINAL-structure indices of the atoms
that survive into the cropped structure, in the order they appear there,
so cropped atom `i` came from original atom `final_indices[i]`. Lets a
caller carry per-atom data across the crop; `polymer_couplings` uses it
to remap real bonds read from a GROMACS `.top` (which describes the
*uncropped* system) onto the cropped structure, so nothing downstream has
to guess bonds from geometry. The mapping is exact because cropping keeps
atoms in sorted original-index order and *repurposes* a deleted heavy
atom as each capping hydrogen rather than appending a new one — so every
cropped atom, caps included, has an original index. Purely additive:
nothing else reads it, and cropping behaviour is unchanged.

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

- **Real perf bug, fixed 2026-08-13: `identify_chains_to_crop` re-fetched
  `universe.atoms.types` (MDAnalysis rebuilds the FULL per-atom type
  array from its internal topology attributes on every access) inside a
  loop over every deleted atom's connections — profiled at ~20s of a
  ~36s crop on a real 58k-atom system (`polymer_couplings`' D18
  reference), scaling with the number of deleted atoms, not a fixed
  cost. Fixed by hoisting the array access out of the loop to a single
  lookup. Two smaller instances of the same "repeated linear-scan
  membership check" pattern fixed alongside it: `if idx in delete`
  checked against a **list** (now a `set`) in both `crop_chains` and
  `_apply_cropping`, and the final atom selection built a
  space-separated `index i1 i2 i3 ...` selection string for MDAnalysis
  to parse (now `universe.atoms[keep_indices]`, direct array indexing).
  Net effect on the same 58k-atom system: ~36s → ~5.3s. Verified
  behaviour-preserving via a new `tests/test_chain_cropper.py` (hand-
  derived keep/delete/replace sets on a small synthetic molecule,
  passing both before and after the fix) and by re-running
  `polymer_couplings`' own non-slow suite (274 passed, unaffected).

- **Real bug, fixed 2026-08-06: the unit cell was dropped.** Both
  `crop_chains` and `_apply_cropping` build the cropped structure with
  `mda.Merge`, which does not carry the box over, and nothing restored
  it. Every cropped structure therefore came out with **no box at all**
  (`PM6_16mer_n20_cropped.gro` ends in `0.00000 0.00000 0.00000` and
  loads as `dimensions=None`), and it propagated: the cropped trajectory
  and every frame sampled from it were unitless too. Any PBC-aware
  calculation downstream was silently non-periodic — which matters
  directly for finding non-covalent contacts across a periodic box.
  Fixed by setting `new_universe.dimensions = universe.dimensions` at
  both `mda.Merge` sites. The parallel trajectory path additionally
  discarded per-frame boxes entirely (`_process_frame` returns only
  coordinates/elements/names), so the box is now collected per frame
  alongside the coordinates and restored before each write. **Per-frame
  boxes are used for trajectory frames, never the structure file's box** —
  under NPT the cell fluctuates, so reusing the structure's box for every
  frame would be wrong. Verified on a synthetic 3-frame trajectory with a
  deliberately different box each frame, serial and parallel paths alike.
  **Note: cropped files produced before this fix carry no box and should
  be regenerated** if anything PBC-aware is to be done with them.

- **Investigated, did not reproduce:** the previously-flagged
  "`add_TopologyAttr` crash in the serial (`-j 1`) trajectory path"
  (repeated `_apply_cropping()` calls on the same `Universe` across
  frames looked like they'd hit MDAnalysis's "attribute already exists"
  error). Live-tested against a real 3-frame slice of a 49960-atom
  trajectory (`poly_workflow/PM6/PM6_traj.trr`) with `-j 1` — completed
  successfully across repeated attribute-add cycles, no crash. Not a
  real bug as far as this test could tell; leave as-is.
- **Fixed 2026-08-06:** `cli.py::generate_output_paths()` — when `-o`
  pointed at a single file (not a directory) *and* both a topology and a
  trajectory were given, the function assigned `-o`'s value straight to
  the topology output slot regardless of its extension (e.g. `-o
  cropped.trr` silently produced two trajectory-shaped outputs and no
  cropped `.gro` at all — found live-testing the same PM6 trajectory
  above; a cropped trajectory with no matching topology can't be opened
  by anything, VMD included, since atom counts have to match). Fixed:
  when both a topology and trajectory are present, both output paths are
  now always derived from the *input* files' own extensions
  (`<stem><input_suffix>` / `<stem>_traj<input_traj_suffix>`), never
  from whatever extension the user's `-o` value happens to have. The
  single-structure case (`-o` given, no trajectory) is unchanged — still
  honors the user's exact filename, since there's no ambiguity there.
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
