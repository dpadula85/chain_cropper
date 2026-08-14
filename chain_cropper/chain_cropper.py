#!/usr/bin/env python3
"""
Enhanced molecular chain cropping tool with trajectory support and configurable parameters.

This module provides functionality to crop alkyl and ether side chains from molecular
structures, with support for trajectories, configurable chain lengths, and multiple
output formats.
"""

import logging
import numpy as np
import MDAnalysis as mda
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from .topology import build_connectivity, side_chain_atoms, sp2_sp3
from .instrumentation import TimingSummary

log = logging.getLogger(__name__)


class ChainCropper:
    """
    A class for cropping molecular side chains from MDAnalysis Universe objects.
    
    This class provides methods to identify and remove alkyl or ether side chains
    while properly capping broken bonds with hydrogen atoms.
    """
    
    def __init__(self, cap_distance: float = 1.09):
        """
        Initialize the ChainCropper.
        
        Parameters
        ----------
        cap_distance : float, default=1.09
            Distance (in Angstrom) for capping hydrogen atoms.
        """
        self.cap_distance = cap_distance
        self.connectivity = None
        # Number of bonds per atom. Read this instead of counting
        # non-negative entries of a `connectivity` row: the row width now
        # follows the most connected atom in the system, so a full row no
        # longer means "four bonds".
        self.degree = None
        self.heavy_atoms = None
        self.sp3_atoms = None
        self.oxygen_atoms = None
        # Boolean heavy-atom mask, kept alongside `connectivity` so the
        # capping routine can tell a deleted heavy neighbour (which becomes
        # a capping hydrogen) from a deleted hydrogen (which just goes).
        self._is_heavy = None
        # See TrajectoryProcessor.final_indices; populated by _build_cropped.
        self.final_indices = None

    def _calculate_new_position(self, anchor_pos: np.ndarray, 
                              old_pos: np.ndarray, distance: float) -> np.ndarray:
        """
        Calculate new position at specified distance from anchor point.
        
        Parameters
        ----------
        anchor_pos : np.ndarray
            Position of the anchor atom
        old_pos : np.ndarray
            Original position of the atom to be moved
        distance : float
            Desired distance from anchor
            
        Returns
        -------
        np.ndarray
            New position coordinates
        """
        bond_vector = old_pos - anchor_pos
        bond_vector = bond_vector / np.linalg.norm(bond_vector)
        return anchor_pos + bond_vector * distance
    
    def _build_connectivity(self, universe: mda.Universe) -> np.ndarray:
        """
        Build the neighbour matrix from bonds, and record the atom degrees.

        Thin wrapper over `chain_cropper.topology.build_connectivity`, kept
        so `self.connectivity`/`self.degree` are populated together.

        Parameters
        ----------
        universe : mda.Universe
            MDAnalysis universe object

        Returns
        -------
        np.ndarray
            Neighbour matrix with -1 as placeholder for empty valence
        """
        self.connectivity, self.degree = build_connectivity(universe)

        return self.connectivity

    def _identify_atom_types(self, universe: mda.Universe) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Identify different atom types in the universe.

        Parameters
        ----------
        universe : mda.Universe
            MDAnalysis universe object

        Returns
        -------
        Tuple[np.ndarray, np.ndarray, np.ndarray]
            Indices of heavy atoms, sp3 atoms, and oxygen atoms
        """
        types = np.asarray(universe.atoms.types, dtype=str)

        self._is_heavy = types != "H"
        heavy_atoms = np.flatnonzero(self._is_heavy)
        _sp2, sp3_atoms = sp2_sp3(
            universe, connectivity=self.connectivity, degree=self.degree
        )
        oxygen_atoms = np.flatnonzero(types == "O")

        return heavy_atoms, sp3_atoms, oxygen_atoms

    def _cap_and_select(self, coords: np.ndarray, elements: np.ndarray,
                        names: np.ndarray, keep: Sequence[int],
                        delete: Sequence[int], replace: Sequence[int],
                        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[int]]:
        """
        Turn each severed bond into a capping hydrogen.

        Single home for the capping step, shared by `crop_chains`,
        `TrajectoryProcessor._apply_cropping` and
        `TrajectoryProcessor._process_frame`. It used to be copy-pasted into
        all three, which is how the O(n) membership test below came to be
        fixed in two of them and left in the third, and how the
        one-cap-per-atom bug below came to be fixed in none.

        A cap does not append a new atom: it *repurposes* the deleted heavy
        atom's own index, moved in along the bond to `cap_distance`. That
        keeps the original K-D bond valid, which is what lets a caller
        remap per-atom data (and real bonds) across the crop via
        `TrajectoryProcessor.final_indices`.

        **Every** deleted heavy neighbour of a kept atom gets its own cap.
        Capping only the first one left a quaternary side-chain carbon that
        lost two substituents -- 9,9-dialkylfluorene and
        indacenodithiophene being the everyday examples -- one hydrogen
        short, with a dangling valence.

        Parameters
        ----------
        coords, elements, names : np.ndarray
            Per-atom arrays for the ORIGINAL (uncropped) atom set. Copied,
            never modified in place.
        keep, delete, replace : sequence of int
            As returned by `identify_chains_to_crop`.

        Returns
        -------
        coords, elements, names : np.ndarray
            Copies with the capping hydrogens written in, still in the
            original atom numbering.
        keep_sorted : list of int
            Original indices of the surviving atoms, ascending, capping
            hydrogens included.
        """
        coords = coords.copy()
        elements = elements.copy()
        names = names.copy()

        atoms_to_keep = set(int(i) for i in keep)
        # Sets, not the incoming lists: this is a per-atom membership test in
        # a loop, so a linear scan of the delete list makes the whole crop
        # quadratic in system size.
        delete_set = set(int(i) for i in delete)

        is_heavy = self._is_heavy
        used: dict = {}

        for atom_idx in replace:
            atom_idx = int(atom_idx)

            for neighbour in self.connectivity[atom_idx]:
                neighbour = int(neighbour)

                if neighbour < 0 or neighbour not in delete_set:
                    continue

                # A deleted hydrogen is not a severed bond to cap -- it went
                # away because its heavy atom did.
                if is_heavy is not None and not is_heavy[neighbour]:
                    continue

                if neighbour in used:
                    raise ValueError(
                        f"cannot cap: deleted atom {neighbour} sits between kept "
                        f"atoms {used[neighbour]} and {atom_idx}, so both need a "
                        f"capping hydrogen but only one atom index is free to "
                        f"carry one. This happens when a saturated ring is cut "
                        f"mid-ring; try a different max_chain_length."
                    )
                used[neighbour] = atom_idx

                coords[neighbour] = self._calculate_new_position(
                    coords[atom_idx], coords[neighbour], self.cap_distance
                )
                elements[neighbour] = 'H'
                names[neighbour] = 'H'
                atoms_to_keep.add(neighbour)

        return coords, elements, names, sorted(atoms_to_keep)


    def identify_chains_to_crop(self, universe: mda.Universe, chain_type: str = 'alkyl',
                               max_chain_length: int = 1) -> Tuple[List[int], List[int], List[int]]:
        """
        Identify atoms to keep, delete, and replace in chain cropping.
        
        Parameters
        ----------
        universe : mda.Universe
            MDAnalysis universe object
        chain_type : str, default='alkyl'
            Type of chains to crop ('alkyl' or 'ether')
        max_chain_length : int, default=1
            Maximum number of heavy atoms to keep in side chains (0 means remove all)
            
        Returns
        -------
        Tuple[List[int], List[int], List[int]]
            Lists of atom indices to keep, delete, and replace with H
        """
        if max_chain_length < 0:
            raise ValueError(
                f"max_chain_length must be non-negative, got {max_chain_length}"
            )

        n_atoms = len(universe.atoms)

        self._build_connectivity(universe)
        self.heavy_atoms, self.sp3_atoms, self.oxygen_atoms = \
            self._identify_atom_types(universe)

        connectivity = self.connectivity
        is_heavy = self._is_heavy

        is_chain = np.zeros(n_atoms, dtype=bool)
        is_chain[side_chain_atoms(universe, chain_type,
                                 connectivity, self.degree)] = True

        # `connectivity` is -1-padded; `safe` lets a padded slot be used as an
        # index into the per-atom masks without tripping negative indexing,
        # and `valid` masks the result back out.
        valid = connectivity >= 0
        safe = np.where(valid, connectivity, 0)

        # A side chain is entered from the core, so the walk starts at every
        # chain atom bonded to a heavy atom that is not itself chain.
        anchored = (valid & is_heavy[safe] & ~is_chain[safe]).any(axis=1)

        # Breadth-first over the chain subgraph, recording how many chain
        # atoms deep each one sits (the anchored atoms are depth 1). This
        # replaces a single-path walk that followed only ONE branch out of
        # each atom: a branched side chain (2-ethylhexyl, or the two alkyls
        # of a 9,9-dialkylfluorene) left the unwalked branch classified as
        # neither kept nor deleted, so it silently vanished from the output
        # while its hydrogens stayed behind, unbonded.
        depth = np.zeros(n_atoms, dtype=int)
        frontier = np.flatnonzero(is_chain & anchored)
        current_depth = 1

        while frontier.size:
            depth[frontier] = current_depth

            neighbours = connectivity[frontier]
            neighbours = np.unique(neighbours[neighbours >= 0])
            frontier = neighbours[is_chain[neighbours] & (depth[neighbours] == 0)]
            current_depth += 1

        # Depth 0 chain atoms were never reached from any core anchor, i.e.
        # they belong to an all-saturated molecule with no core to hang off
        # (a solvent alkane, say). There is no side chain to trim there, so
        # they are kept -- the old code dropped them.
        delete_mask = is_chain & (depth > max_chain_length)

        # Hydrogens of a deleted heavy atom go with it. Reading the heavy
        # mask off `self._is_heavy` rather than `universe.atoms.types` keeps
        # this off MDAnalysis's per-access attribute rebuild, which used to
        # dominate the runtime of a large crop.
        deleted = np.flatnonzero(delete_mask)
        if deleted.size:
            neighbours = connectivity[deleted]
            neighbours = np.unique(neighbours[neighbours >= 0])
            delete_mask[neighbours[~is_heavy[neighbours]]] = True

        keep_mask = ~delete_mask

        # A kept atom needs one capping hydrogen per deleted heavy neighbour.
        # The old code only ever considered kept *chain* atoms here, so
        # `max_chain_length=0` -- documented as "remove all side chains" --
        # deleted every chain atom and then capped nothing at all, handing
        # back a core full of dangling valences.
        deleted_heavy = delete_mask & is_heavy
        needs_cap = (valid & deleted_heavy[safe]).any(axis=1)
        replace_atoms = np.flatnonzero(keep_mask & needs_cap)

        # Each cap is carried by the deleted atom's own index, so a deleted
        # atom bridging two kept atoms cannot cap both. Fail loudly rather
        # than emit a structure that is quietly one hydrogen short.
        bridging = np.flatnonzero(
            deleted_heavy & ((valid & keep_mask[safe]).sum(axis=1) > 1)
        )
        if bridging.size:
            raise ValueError(
                f"cannot cap: {bridging.size} deleted atom(s) sit between two "
                f"kept atoms (first is atom {int(bridging[0])}), so each would "
                f"have to become two capping hydrogens. This happens when a "
                f"saturated ring is cut mid-ring; try a different "
                f"max_chain_length ({max_chain_length} was requested)."
            )

        return (np.flatnonzero(keep_mask).tolist(),
                np.flatnonzero(delete_mask).tolist(),
                replace_atoms.tolist())

    def crop_chains(self, universe: mda.Universe, chain_type: str = 'alkyl',
                   max_chain_length: int = 1) -> Tuple[mda.Universe, List[int], List[int]]:
        """
        Crop side chains from a universe.
        
        Parameters
        ----------
        universe : mda.Universe
            Input universe
        chain_type : str, default='alkyl'
            Type of chains to crop
        max_chain_length : int, default=1
            Maximum chain length to keep

        Returns
        -------
        Tuple[mda.Universe, List[int], List[int]]
            New universe with cropped chains, followed by the `keep` and
            `replace` atom-index lists from `identify_chains_to_crop`
            (the atoms retained, and those capped with a new H).

            The input universe is left untouched. It used to have its
            positions, types and names overwritten in place, so a caller
            that cropped and then went back to the original got the cropped
            atom types instead.
        """
        keep, delete, replace = self.identify_chains_to_crop(
            universe, chain_type, max_chain_length
        )

        return self._build_cropped(universe, keep, delete, replace), keep, replace

    def _build_cropped(self, universe: mda.Universe, keep: Sequence[int],
                       delete: Sequence[int], replace: Sequence[int],
                       ) -> mda.Universe:
        """
        Assemble the cropped universe: cap the severed bonds, keep the
        surviving atoms, carry the unit cell over.

        Shared by `crop_chains` and `TrajectoryProcessor._apply_cropping` so
        the two cannot drift apart -- they previously wrote different
        topology attributes, `crop_chains` setting `types` and
        `_apply_cropping` setting `name`/`element`, so only one of the two
        returned a universe with a usable `atoms.elements`.

        Parameters
        ----------
        universe : mda.Universe
            Source universe, read-only. Its currently loaded frame supplies
            the coordinates and the box.
        keep, delete, replace : sequence of int
            As returned by `identify_chains_to_crop`.

        Returns
        -------
        mda.Universe
            The cropped structure, with `name`, `type` and `element` set.
        """
        coords, elements, names, keep_indices = self._cap_and_select(
            universe.atoms.positions,
            np.asarray(universe.atoms.types, dtype=object),
            np.asarray(universe.atoms.names, dtype=object),
            keep, delete, replace,
        )

        # Record the cropped -> original index mapping (see
        # TrajectoryProcessor.final_indices).
        self.final_indices = np.asarray(keep_indices)

        # Direct array indexing -- avoids building and parsing a selection
        # string with one token per atom, which does not scale to large
        # systems (tens of thousands of atoms).
        new_universe = mda.Merge(universe.atoms[keep_indices])

        # Write the capped positions/labels onto the OUTPUT universe rather
        # than onto the input, so cropping has no side effect on its
        # argument.
        new_universe.atoms.positions = coords[keep_indices]
        kept_names = np.asarray([names[i] for i in keep_indices], dtype=str)
        kept_elements = np.asarray([elements[i] for i in keep_indices], dtype=str)
        new_universe.add_TopologyAttr('name', kept_names)
        new_universe.add_TopologyAttr('type', kept_elements)
        new_universe.add_TopologyAttr('element', kept_elements)

        # mda.Merge does not carry the unit cell over, so without this the
        # cropped structure comes out with no box at all and every
        # PBC-aware calculation downstream silently becomes non-periodic.
        # universe.dimensions is the box of the frame currently loaded, so
        # this picks up that frame's box when called from the trajectory
        # loop, and the structure's own box when called on a structure.
        new_universe.dimensions = universe.dimensions

        return new_universe


class TrajectoryProcessor(ChainCropper):
    """
    Handler for processing trajectories with chain cropping.
    Inherits from ChainCropper and reuses cropping indices across frames.
    Supports parallel processing for trajectory files.
    """
    
    def __init__(self, cap_distance: float = 1.09):
        """
        Initialize the TrajectoryProcessor.
        
        Parameters
        ----------
        cap_distance : float, default=1.09
            Distance (in Angstrom) for capping hydrogen atoms.
        """
        super().__init__(cap_distance)
        self.keep_indices = None
        self.replace_indices = None
        self.delete_indices = None
        # Indices, in the ORIGINAL structure's numbering, of the atoms that
        # survive into the cropped structure, in the order they appear there.
        # So cropped atom i corresponds to original atom final_indices[i].
        # Populated by _apply_cropping; lets a caller carry per-atom data
        # (e.g. real bonds read from a GROMACS topology) across the crop.
        self.final_indices = None
        # Identity of the analysis currently cached in
        # keep/delete/replace_indices -- see `analyze`.
        self._analysis_params = None
        self._analyzed_structure_path = None
        self._analyzed_universe = None
        self.timing = TimingSummary()
    
    def analyze(self, structure_file: Optional[str] = None,
                chain_type: str = 'alkyl', max_chain_length: int = 1,
                structure_universe: Optional[mda.Universe] = None) -> mda.Universe:
        """
        Determine which atoms to keep/delete/cap, caching the result.

        A repeat call with the same structure identity, chain_type and
        max_chain_length returns the cached universe and skips
        `identify_chains_to_crop` -- the expensive step, dominated by bond
        guessing on a large system. Keyed on all three, not just "has this
        run once", so reusing one processor with a different
        `max_chain_length` still re-analyses instead of serving stale
        indices.

        Parameters
        ----------
        structure_file : Optional[str]
            Input structure file. Ignored if `structure_universe` is given.
        chain_type : str, default='alkyl'
        max_chain_length : int, default=1
        structure_universe : Optional[MDAnalysis.Universe], default=None
            Pass one carrying real bonds instead of loading `structure_file`
            -- see `process_trajectory`'s docstring for why.

        Returns
        -------
        mda.Universe
            The universe the analysis ran against.
        """
        if structure_universe is None and structure_file is None:
            raise ValueError("Must supply structure_file or structure_universe")

        params = (chain_type, max_chain_length)
        resolved_path = (
            str(Path(structure_file).resolve()) if structure_file is not None else None
        )

        # A hit needs the SAME universe object previously analyzed, identified
        # by id (if handed one directly) or by its load path (if not) -- a
        # caller can reasonably do either on a repeat call, e.g. cli.py loads
        # by path once, then passes the returned object back in explicitly.
        cached = (self._analyzed_universe is not None and params == self._analysis_params)
        if cached and structure_universe is not None:
            cached = id(structure_universe) == id(self._analyzed_universe)
        elif cached:
            cached = resolved_path == self._analyzed_structure_path

        if cached:
            log.debug("Reusing cached analysis for %s", params)
            return self._analyzed_universe

        if structure_universe is None:
            with self.timing.measure("structure_load"):
                log.info("Loading structure from %s...", structure_file)
                structure_universe = mda.Universe(structure_file)
        else:
            log.info("Using the supplied structure universe for connectivity")

        with self.timing.measure("analysis"):
            log.info("Analyzing chain structure...")
            keep, delete, replace = self.identify_chains_to_crop(
                structure_universe, chain_type, max_chain_length
            )

        self.keep_indices = keep
        self.delete_indices = delete
        self.replace_indices = replace
        self._analysis_params = params
        self._analyzed_structure_path = resolved_path
        self._analyzed_universe = structure_universe

        log.info("Identified %d atoms to keep, %d to delete, %d to cap",
                  len(keep), len(delete), len(replace))

        return structure_universe

    def write_structure(self, universe: mda.Universe, output_path: Path) -> mda.Universe:
        """
        Crop `universe` with the currently cached indices and write it out.

        Parameters
        ----------
        universe : mda.Universe
            Universe to crop (must have the same atom count as the one
            `analyze` ran against).
        output_path : Path
            Where to write the cropped structure.

        Returns
        -------
        mda.Universe
            The cropped structure.
        """
        with self.timing.measure("structure_write"):
            cropped = self._apply_cropping(universe)
            log.info("Writing cropped structure to %s...", output_path)
            cropped.atoms.write(str(output_path))

        return cropped

    def process_trajectory(self, structure_file: str, output_path: str,
                          trajectory_file: Optional[str] = None,
                          chain_type: str = 'alkyl',
                          max_chain_length: int = 1,
                          n_jobs: int = -1,
                          structure_universe: Optional[mda.Universe] = None) -> None:
        """
        Process trajectory and write to file.

        This method crops the structure once to determine which atoms to keep,
        then applies the same cropping to all frames in the trajectory.

        Parameters
        ----------
        structure_file : str
            Input structure file (e.g., .gro, .pdb)
        output_path : str
            Output file path
        trajectory_file : Optional[str], default=None
            Trajectory file (e.g., .xtc, .trr, .dcd). If None, processes only structure.
        chain_type : str, default='alkyl'
            Type of chains to crop ('alkyl' or 'ether')
        max_chain_length : int, default=1
            Maximum chain length to keep
        n_jobs : int, default=-1
            Number of parallel jobs for trajectory processing.
            -1 uses all available cores, 1 disables parallelization
        structure_universe : Optional[MDAnalysis.Universe], default=None
            Universe to take the connectivity from, instead of loading
            `structure_file`. Pass one carrying real bonds -- read from a
            GROMACS topology, say -- so which atoms count as side chain is
            decided from the true connectivity rather than from bonds
            guessed off the geometry. Guessing can both miss real bonds and
            invent absent ones on a strained snapshot, and it is also the
            slow step on a large system. `structure_file` is still used for
            the coordinates and for pairing with the trajectory, so the two
            must describe the same atoms in the same order.
        """
        structure_universe = self.analyze(
            structure_file=structure_file, chain_type=chain_type,
            max_chain_length=max_chain_length, structure_universe=structure_universe,
        )

        output_path = Path(output_path)

        # If no trajectory, just write the structure
        if trajectory_file is None:
            self.write_structure(structure_universe, output_path)
            log.info("Done!")
            return

        # Load universe with trajectory
        with self.timing.measure("trajectory_load"):
            log.info("Loading trajectory from %s...", trajectory_file)
            traj_universe = mda.Universe(structure_file, trajectory_file)
            n_frames = len(traj_universe.trajectory)
        log.info("Processing %d frames...", n_frames)

        extension = output_path.suffix.lower()

        # Write trajectory or single frame
        if extension in ['.xtc', '.trr', '.dcd']:
            self._write_trajectory(traj_universe, output_path, n_frames, n_jobs)
        elif extension in ['.xyz', '.gro', '.pdb']:
            # Just write the first frame
            traj_universe.trajectory[0]
            self.write_structure(traj_universe, output_path)
            log.info("Done!")
        else:
            raise ValueError(f"Unsupported output format: {extension}")

    def _apply_cropping(self, universe: mda.Universe) -> mda.Universe:
        """
        Apply stored cropping indices to a universe.
        
        Parameters
        ----------
        universe : mda.Universe
            Universe to crop (must have same atom count as original)
            
        Returns
        -------
        mda.Universe
            Cropped universe with capped bonds
        """
        if self.keep_indices is None:
            raise RuntimeError("Must run identify_chains_to_crop first")

        return self._build_cropped(
            universe, self.keep_indices, self.delete_indices, self.replace_indices
        )


    def _process_frame(self, frame_idx: int, coords: np.ndarray, 
                      elements: np.ndarray, names: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Process a single frame's coordinates.
        
        Parameters
        ----------
        frame_idx : int
            Frame index (for logging/debugging)
        coords : np.ndarray
            Coordinates for this frame
        elements : np.ndarray
            Element types for all atoms
        names : np.ndarray
            Atom names for all atoms
            
        Returns
        -------
        Tuple[np.ndarray, np.ndarray, np.ndarray]
            Cropped (coordinates, elements, names) for atoms to keep
        """
        coords, elements, names, keep_indices = self._cap_and_select(
            coords, elements, names,
            self.keep_indices, self.delete_indices, self.replace_indices,
        )

        return (coords[keep_indices],
                elements[keep_indices],
                names[keep_indices])
    
    def _write_trajectory(self, universe: mda.Universe, output_path: Path,
                         n_frames: int, n_jobs: int = -1) -> None:
        """
        Write trajectory file with parallel processing and progress bar.
        
        Parameters
        ----------
        universe : mda.Universe
            Universe with trajectory
        output_path : Path
            Output file path
        n_frames : int
            Number of frames to process
        n_jobs : int, default=-1
            Number of parallel jobs. -1 uses all available cores,
            1 disables parallelization
        """
        from tqdm import tqdm
        from tqdm.contrib.logging import logging_redirect_tqdm
        from joblib import Parallel, delayed
        import multiprocessing as mp

        # Process first frame to get structure info
        universe.trajectory[0]
        first_cropped = self._apply_cropping(universe)
        n_atoms = len(first_cropped.atoms)

        # Determine number of jobs
        if n_jobs == -1:
            n_jobs = mp.cpu_count()
        elif n_jobs < 1:
            n_jobs = 1

        log.info("Processing %d frames using %d core(s)...", n_frames, n_jobs)

        if n_jobs == 1:
            # Serial processing
            with mda.Writer(str(output_path), n_atoms) as writer, \
                    logging_redirect_tqdm():
                for ts in tqdm(universe.trajectory, total=n_frames,
                              desc="Writing trajectory", unit="frame"):
                    with self.timing.measure("per_frame_cropping"):
                        cropped_frame = self._apply_cropping(universe)
                    with self.timing.measure("trajectory_write"):
                        writer.write(cropped_frame.atoms)
        else:
            # Parallel processing
            # First, load all coordinates into memory
            log.info("Loading trajectory into memory...")
            all_coords = []
            all_dimensions = []
            elements = universe.atoms.types.copy()
            names = universe.atoms.names.copy()

            with self.timing.measure("trajectory_load"), logging_redirect_tqdm():
                for ts in tqdm(universe.trajectory, total=n_frames,
                              desc="Loading frames", unit="frame"):
                    all_coords.append(universe.atoms.positions.copy())
                    # Each frame keeps its own box: under NPT the cell
                    # fluctuates, so the structure file's box must not be
                    # reused for every frame.
                    all_dimensions.append(
                        None if ts.dimensions is None else ts.dimensions.copy()
                    )

            # Process frames in parallel
            log.info("Processing frames in parallel with %d workers...", n_jobs)
            with self.timing.measure("per_frame_cropping"), logging_redirect_tqdm():
                processed_frames = Parallel(n_jobs=n_jobs)(
                    delayed(self._process_frame)(i, coords, elements, names)
                    for i, coords in enumerate(tqdm(all_coords,
                                                   desc="Processing frames",
                                                   unit="frame"))
                )

            # Write processed frames
            log.info("Writing trajectory...")
            with self.timing.measure("trajectory_write"), \
                    mda.Writer(str(output_path), n_atoms) as writer, \
                    logging_redirect_tqdm():
                # Create temporary universe for writing
                temp_u = mda.Universe.empty(n_atoms,
                                           n_residues=1,
                                           atom_resindex=[0]*n_atoms,
                                           trajectory=True)
                temp_u.add_TopologyAttr('type', processed_frames[0][1])
                temp_u.add_TopologyAttr('name', processed_frames[0][2])

                for (coords, elements, names), dimensions in tqdm(
                        zip(processed_frames, all_dimensions),
                        total=len(processed_frames),
                        desc="Writing frames", unit="frame"):
                    temp_u.atoms.positions = coords
                    temp_u.atoms.types = elements
                    temp_u.atoms.names = names
                    # Restore this frame's own box before writing it.
                    if dimensions is not None:
                        temp_u.dimensions = dimensions
                    writer.write(temp_u.atoms)

        log.info("Trajectory written to %s", output_path)


# Example usage
if __name__ == '__main__':
    # Create processor
    processor = TrajectoryProcessor(cap_distance=1.09)
    
    # Process structure only
    processor.process_trajectory(
        structure_file='system.gro',
        output_path='cropped.gro',
        chain_type='alkyl',
        max_chain_length=1
    )
    
    # Process structure + trajectory with all CPU cores
    processor = TrajectoryProcessor(cap_distance=1.09)
    processor.process_trajectory(
        structure_file='system.gro',
        trajectory_file='traj.xtc',
        output_path='cropped_traj.xtc',
        chain_type='alkyl',
        max_chain_length=1,
        n_jobs=-1  # Use all cores
    )
    
    # Process with 4 cores
    processor = TrajectoryProcessor(cap_distance=1.09)
    processor.process_trajectory(
        structure_file='system.gro',
        trajectory_file='traj.xtc',
        output_path='cropped_traj.xtc',
        chain_type='alkyl',
        max_chain_length=1,
        n_jobs=4
    )
    
    # Serial processing (no parallelization)
    processor = TrajectoryProcessor(cap_distance=1.09)
    processor.process_trajectory(
        structure_file='system.gro',
        trajectory_file='traj.xtc',
        output_path='cropped_traj.xtc',
        chain_type='alkyl',
        max_chain_length=1,
        n_jobs=1
    )
