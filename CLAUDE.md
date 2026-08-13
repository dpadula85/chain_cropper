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

## Connectivity/hybridisation module

`chain_cropper/topology.py` is the single home for the bond ->
neighbour-matrix -> sp2/sp3 classification chain: `build_connectivity`,
`sp2_sp3`/`get_sp2`, `side_chain_atoms`. `SelIntCoords.sel_intcoords` and
`oligomer_builder.enhanced_breaker` both import `get_sp2` from here rather
than carrying their own copy — see "Fixed 2026-08-13" below for why that
mattered. `oligomer_builder` and `SelIntCoords` both now depend on this
package (`chain-cropper` in their `install_requires`).

`side_chain_atoms` (and therefore `identify_chains_to_crop`) treats a
saturated atom as side chain only if it has **at most one** core (non-
side-chain heavy) neighbour. An atom with two or more is a *bridge* —
holding two pieces of core together — and counts as core itself. This is
what keeps fluorene/cyclopentadithiophene/indacenodithiophene-type fused
units intact: their sp3 bridge carbon is four-coordinate, same as a plain
alkyl chain carbon, so without this rule it reads as side chain too.

## Known gaps

- No CI beyond the GitHub Actions workflow added 2026-08-13 (pytest
  across Python 3.9-3.12); no live cluster smoke-test in CI, so a real
  MD trajectory still needs a manual run to catch anything the synthetic
  test molecules don't cover.

## Fixed 2026-08-13

Six real defects, all confirmed against both synthetic molecules and
real production `.gro` files (PM6/D18/DPP/pg2TT from `poly_workflow`)
before and after:

1. **Fused sp3 bridge carbons silently vanished or split their ring
   system.** `identify_chains_to_crop` treated every four-coordinate atom
   as side chain, including fluorene/CPDT/IDT-type bridge carbons that
   hold two ring systems together. Confirmed on real `pg2TT` data: 96
   heavy atoms (its ether-linked terminal methyls, a related but
   distinct case — see #2) were dropped from a single frame. Fixed by
   `side_chain_atoms`' bridge-vs-pendant rule (see above).
2. **`max_chain_length=0` ("remove all side chains") capped nothing.**
   The capping loop only ever considered kept *chain* atoms; at
   `max_chain_length=0` there are none, so every anchor was left with a
   dangling valence and no capping hydrogen. Fixed: the new capping step
   caps every deleted heavy neighbour of every *kept* atom, chain or not.
3. **A branched side chain lost its unwalked branch.** The single-path
   walk (`_find_chain_from_terminus`) followed exactly one neighbour out
   of each chain atom; a second branch (2-ethylhexyl, a
   9,9-dialkylfluorene's second arm) was left classified as neither kept
   nor deleted — its heavy atom vanished from the output while its
   hydrogens stayed behind, bonded to nothing. Fixed by a breadth-first,
   depth-limited walk (`identify_chains_to_crop`) that reaches every
   branch.
4. **A quaternary/gem-disubstituted atom got only one capping hydrogen.**
   The capping loop used `deleted_connected[0]` — the first severed bond
   only. A bridge or side-chain carbon losing two substituents at once
   came out one hydrogen short, silently. Fixed: `_cap_and_select` caps
   every severed heavy bond of a kept atom, and raises rather than
   silently dropping a cap if two kept atoms would need to share one
   deleted atom's index (a saturated ring cut mid-ring).
5. **`_build_connectivity` hard-coded 4 bonds/atom and dropped the
   5th+.** Fixed: `topology.build_connectivity` widens to the true
   maximum degree.
6. **`crop_chains` mutated its input Universe and returned a
   differently-shaped one than `_apply_cropping`.** It overwrote the
   input's positions/types/names in place, and set `type` where
   `_apply_cropping` set `name`/`element` (so only one of the two callers
   got a working `atoms.elements`). The parallel trajectory path
   (`_process_frame`) also still did an O(n) list-membership scan already
   fixed in the serial path. Fixed: both entry points and the trajectory
   frame path now share one capping-and-selection routine
   (`_cap_and_select`/`_build_cropped`), which copies rather than
   mutates its input.

`ether=True`'s no-op branch (previously listed here as a known gap) is
fixed as a side effect of #1's rewrite — see the topology module note
above; `SelIntCoords`/`oligomer_builder`'s copies got the same fix by
switching to import this package's implementation instead of carrying
their own.

## Dependencies

MDAnalysis>=2.0.0, numpy>=1.19.0, tqdm, joblib. No vendored code.
