#!/usr/bin/env python3
"""
Enhanced molecular chain cropping tool with trajectory support and configurable parameters.

This module provides functionality to crop alkyl and ether side chains from molecular
structures, with support for trajectories, configurable chain lengths, and multiple
output formats.
"""

import argparse
import numpy as np
import MDAnalysis as mda
from pathlib import Path
from typing import List, Optional, Tuple, Union, Dict, Set


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
        self.heavy_atoms = None
        self.sp3_atoms = None
        self.oxygen_atoms = None
    
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
        Build connectivity matrix from bonds.
        
        Parameters
        ----------
        universe : mda.Universe
            MDAnalysis universe object
            
        Returns
        -------
        np.ndarray
            Connectivity matrix with -1 as placeholder for empty valence
        """
        # Get or guess bonds
        try:
            bonds = universe.bonds.to_indices()
        except (AttributeError, ValueError):
            universe.atoms.guess_bonds()
            bonds = universe.bonds.to_indices()
        
        # Initialize connectivity matrix
        max_bonds = 4  # Assume maximum 4 bonds per atom
        connectivity = np.full((len(universe.atoms), max_bonds), -1, dtype=int)
        
        # Fill connectivity matrix
        for bond in bonds:
            at1, at2 = bond
            for j in range(max_bonds):
                if connectivity[at1, j] == -1:
                    connectivity[at1, j] = at2
                    break
            
            for j in range(max_bonds):
                if connectivity[at2, j] == -1:
                    connectivity[at2, j] = at1
                    break
        
        return connectivity
    
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
        heavy_atoms = np.where(universe.atoms.types != "H")[0]
        sp3_atoms = np.where(np.all(self.connectivity > -1, axis=1))[0]
        oxygen_atoms = np.where(universe.atoms.types == "O")[0]
        
        return heavy_atoms, sp3_atoms, oxygen_atoms
    
    def _find_chain_from_terminus(self, start_atom: int, chain_type: str, 
                                 max_length: int) -> List[int]:
        """
        Find a chain starting from a terminal atom.
        
        Parameters
        ----------
        start_atom : int
            Starting atom index
        chain_type : str
            Type of chain to find
        max_length : int
            Maximum chain length
            
        Returns
        -------
        List[int]
            List of atom indices in the chain
        """
        if chain_type == 'alkyl':
            valid_atoms = set(self.sp3_atoms)
        elif chain_type == 'ether':
            valid_atoms = set(np.concatenate([self.sp3_atoms, self.oxygen_atoms]))
        else:
            valid_atoms = set(self.sp3_atoms)
        
        chain = [start_atom]
        current = start_atom
        visited = {start_atom}
        
        while True:
            # Find next atom in chain
            connected = self.connectivity[current]
            connected_valid = [idx for idx in connected 
                             if idx >= 0 and idx in valid_atoms and idx not in visited]
           
            if not connected_valid:
                break
            
            # Choose the first valid connected atom
            next_atom = connected_valid[0]

            if next_atom in visited:
                break 

            chain.append(next_atom)
            visited.add(next_atom)
            current = next_atom
        
        return chain
    
    def _recursive_delete(self, atom_idx: int, keep_atoms: Set[int], 
                         delete_atoms: Set[int], visited: Set[int]) -> None:
        """
        Recursively find all atoms connected to atom_idx that should be deleted.
        
        Parameters
        ----------
        atom_idx : int
            Current atom index
        keep_atoms : Set[int]
            Set of atoms that must be kept
        delete_atoms : Set[int]
            Set of atoms to delete (modified in place)
        visited : Set[int]
            Set of already visited atoms (modified in place)
        """
        if atom_idx in visited or atom_idx in keep_atoms:
            return
            
        visited.add(atom_idx)
        delete_atoms.add(atom_idx)
        
        # Get all connected atoms
        connected = self.connectivity[atom_idx]
        for connected_atom in connected:
            if connected_atom >= 0:  # Valid connection (not -1)
                self._recursive_delete(connected_atom, keep_atoms, delete_atoms, visited)
    
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
        self.connectivity = self._build_connectivity(universe)
        self.heavy_atoms, self.sp3_atoms, self.oxygen_atoms = self._identify_atom_types(universe)
        
        if chain_type == 'alkyl':
            chain_atoms = self.sp3_atoms
        elif chain_type == 'ether':
            chain_atoms = np.concatenate([self.sp3_atoms, self.oxygen_atoms])
        else:
            chain_atoms = self.sp3_atoms
        
        keep_chain_atoms = []
        delete_chain_atoms = []
        
        # Find terminal atoms (connected to non-chain atoms)
        terminal_atoms = []
        for atom_idx in chain_atoms:
            connected = self.connectivity[atom_idx]
            connected_heavy = connected[np.isin(connected, self.heavy_atoms)]
            connected_heavy = connected_heavy[connected_heavy >= 0]
            
            # Check if connected to non-chain atoms
            non_chain_connected = [idx for idx in connected_heavy if idx not in chain_atoms]
            if non_chain_connected:
                terminal_atoms.append(atom_idx)
        
        # Process each terminal atom to identify chain atoms to keep
        processed = set()
        for terminal in terminal_atoms:
            if terminal in processed:
                continue
                
            # Find chain from this terminal
            chain = self._find_chain_from_terminus(terminal, chain_type, max_chain_length + 1)
            
            if len(chain) > max_chain_length:
                # Keep first max_chain_length atoms
                keep_chain_atoms.extend(chain[:max_chain_length])
                delete_chain_atoms.extend(chain[max_chain_length:])
            
            processed.update(chain)

        # All atoms that should be kept (core structure + chain atoms to keep)
        all_atoms = set(range(len(universe.atoms)))
        keep_atoms = set(keep_chain_atoms)
        
        # Add all non-chain atoms to keep set
        for atom_idx in all_atoms:
            if atom_idx not in chain_atoms:
                keep_atoms.add(atom_idx)
        
        # Find atoms that need to be replaced (atoms in keep set that are connected to atoms not in keep set)
        replace_atoms = []
        delete_atoms = set()
        
        # Start with the identified chain atoms to delete
        atoms_to_delete = set(delete_chain_atoms)
        
        # Recursively find all atoms connected to the delete_chain_atoms that should also be deleted
        for atom_idx in delete_chain_atoms:
            visited = set()
            self._recursive_delete(atom_idx, keep_atoms, atoms_to_delete, visited)
        
        # Find replacement points - atoms in keep set connected to deleted atoms
        for atom_idx in keep_chain_atoms:
            connected = self.connectivity[atom_idx]
            connected_atoms = [idx for idx in connected if idx >= 0]
            
            # Check if connected to any atom that will be deleted
            has_deletable_connection = any(conn in atoms_to_delete for conn in connected_atoms)
            
            if has_deletable_connection:
                replace_atoms.append(atom_idx)
        
        # Add all hydrogens connected to deleted heavy atoms
        all_atoms_to_delete = atoms_to_delete.copy()
        for atom_idx in atoms_to_delete:
            connected = self.connectivity[atom_idx]
            for connected_atom in connected:
                if connected_atom >= 0 and universe.atoms.types[connected_atom] == "H":
                    all_atoms_to_delete.add(connected_atom)
        
        # Final keep set excludes all deleted atoms
        final_keep_atoms = keep_atoms - all_atoms_to_delete
        
        return list(final_keep_atoms), list(all_atoms_to_delete), replace_atoms

    def crop_chains(self, universe: mda.Universe, chain_type: str = 'alkyl',
                   max_chain_length: int = 1) -> mda.Universe:
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
        mda.Universe
            New universe with cropped chains
        """
        keep, delete, replace = self.identify_chains_to_crop(
            universe, chain_type, max_chain_length
        )
        
        # Determine which atoms to keep in final structure
        atoms_to_keep = set(keep)
        
        # Cap broken bonds
        coords = universe.atoms.positions.copy()
        elements = universe.atoms.types.copy()
        names = universe.atoms.names.copy()
        
        for atom_idx in replace:
            # Find the deleted atom it was connected to
            connected = self.connectivity[atom_idx]
            connected_atoms = [idx for idx in connected if idx >= 0]
            
            deleted_connected = [idx for idx in connected_atoms if idx in delete]
            if deleted_connected:
                # Use the first deleted connection for positioning
                deleted_atom = deleted_connected[0]
                anchor_pos = coords[atom_idx]
                old_pos = coords[deleted_atom]
                
                # Create new hydrogen at the deleted atom position
                new_h_pos = self._calculate_new_position(anchor_pos, old_pos, self.cap_distance)
                coords[deleted_atom] = new_h_pos
                elements[deleted_atom] = 'H'
                names[deleted_atom] = 'H'
                
                # Add the hydrogen back to atoms to keep
                atoms_to_keep.add(deleted_atom)
        
        # Create selection string
        keep_indices = sorted(list(atoms_to_keep))
        
        # Update universe
        universe.atoms.positions = coords
        universe.atoms.types = elements
        universe.atoms.names = names
        
        # Create new universe with selected atoms
        selection = universe.select_atoms(f'index {" ".join(map(str, keep_indices))}')
        new_universe = mda.Merge(selection)
        
        return new_universe, keep, replace


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
    
    def process_trajectory(self, structure_file: str, output_path: str,
                          trajectory_file: Optional[str] = None,
                          chain_type: str = 'alkyl', 
                          max_chain_length: int = 1,
                          n_jobs: int = -1) -> None:
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
        """
        from tqdm import tqdm
        
        # Load structure to determine cropping indices
        print(f"Loading structure from {structure_file}...")
        structure_universe = mda.Universe(structure_file)
        
        # Perform cropping analysis on structure only
        print("Analyzing chain structure...")
        keep, delete, replace = self.identify_chains_to_crop(
            structure_universe, chain_type, max_chain_length
        )
        
        # Store indices for reuse
        self.keep_indices = keep
        self.delete_indices = delete
        self.replace_indices = replace
        
        print(f"Identified {len(keep)} atoms to keep, {len(delete)} to delete, {len(replace)} to cap")
        
        # Create template cropped structure
        cropped_structure = self._apply_cropping(structure_universe)
        
        # Determine output format
        output_path = Path(output_path)
        extension = output_path.suffix.lower()
        
        # If no trajectory, just write the structure
        if trajectory_file is None:
            print(f"Writing cropped structure to {output_path}...")
            cropped_structure.atoms.write(str(output_path))
            print("Done!")
            return
        
        # Load universe with trajectory
        print(f"Loading trajectory from {trajectory_file}...")
        traj_universe = mda.Universe(structure_file, trajectory_file)
        n_frames = len(traj_universe.trajectory)
        print(f"Processing {n_frames} frames...")
        
        # Write trajectory or single frame
        if extension in ['.xtc', '.trr', '.dcd']:
            self._write_trajectory(traj_universe, output_path, n_frames, n_jobs)
        elif extension in ['.xyz', '.gro', '.pdb']:
            # Just write the first frame
            traj_universe.trajectory[0]
            cropped_frame = self._apply_cropping(traj_universe)
            cropped_frame.atoms.write(str(output_path))
            print("Done!")
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
        
        # Work with copies
        coords = universe.atoms.positions.copy()
        elements = universe.atoms.types.copy()
        names = universe.atoms.names.copy()
        
        atoms_to_keep = set(self.keep_indices)
        
        # Cap broken bonds
        for atom_idx in self.replace_indices:
            # Find the deleted atom it was connected to
            connected = self.connectivity[atom_idx]
            connected_atoms = [idx for idx in connected if idx >= 0]
            
            deleted_connected = [idx for idx in connected_atoms 
                               if idx in self.delete_indices]
            if deleted_connected:
                # Use the first deleted connection for positioning
                deleted_atom = deleted_connected[0]
                anchor_pos = coords[atom_idx]
                old_pos = coords[deleted_atom]
                
                # Create new hydrogen at the deleted atom position
                new_h_pos = self._calculate_new_position(
                    anchor_pos, old_pos, self.cap_distance
                )
                coords[deleted_atom] = new_h_pos
                elements[deleted_atom] = 'H'
                names[deleted_atom] = 'H'
                
                # Add the hydrogen back to atoms to keep
                atoms_to_keep.add(deleted_atom)
        
        # Create selection
        keep_indices_sorted = sorted(list(atoms_to_keep))
        
        # Update universe temporarily
        universe.atoms.positions = coords
        universe.atoms.types = elements
        universe.atoms.names = names
        
        # Create new universe with selected atoms
        selection = universe.select_atoms(
            f'index {" ".join(map(str, keep_indices_sorted))}'
        )
        new_universe = mda.Merge(selection)
        
        return new_universe
    
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
        atoms_to_keep = set(self.keep_indices)
        coords_copy = coords.copy()
        elements_copy = elements.copy()
        names_copy = names.copy()
        
        # Cap broken bonds
        for atom_idx in self.replace_indices:
            connected = self.connectivity[atom_idx]
            connected_atoms = [idx for idx in connected if idx >= 0]
            
            deleted_connected = [idx for idx in connected_atoms 
                               if idx in self.delete_indices]
            if deleted_connected:
                deleted_atom = deleted_connected[0]
                anchor_pos = coords_copy[atom_idx]
                old_pos = coords_copy[deleted_atom]
                
                # Create new hydrogen position
                new_h_pos = self._calculate_new_position(
                    anchor_pos, old_pos, self.cap_distance
                )
                coords_copy[deleted_atom] = new_h_pos
                elements_copy[deleted_atom] = 'H'
                names_copy[deleted_atom] = 'H'
                atoms_to_keep.add(deleted_atom)
        
        # Return only data for atoms to keep
        keep_indices_sorted = sorted(list(atoms_to_keep))
        return (coords_copy[keep_indices_sorted], 
                elements_copy[keep_indices_sorted],
                names_copy[keep_indices_sorted])
    
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
        
        print(f"Processing {n_frames} frames using {n_jobs} core(s)...")
        
        if n_jobs == 1:
            # Serial processing
            with mda.Writer(str(output_path), n_atoms) as writer:
                for ts in tqdm(universe.trajectory, total=n_frames, 
                              desc="Writing trajectory", unit="frame"):
                    cropped_frame = self._apply_cropping(universe)
                    writer.write(cropped_frame.atoms)
        else:
            # Parallel processing
            # First, load all coordinates into memory
            print("Loading trajectory into memory...")
            all_coords = []
            elements = universe.atoms.types.copy()
            names = universe.atoms.names.copy()
            
            for ts in tqdm(universe.trajectory, total=n_frames, 
                          desc="Loading frames", unit="frame"):
                all_coords.append(universe.atoms.positions.copy())
            
            # Process frames in parallel
            print(f"Processing frames in parallel with {n_jobs} workers...")
            processed_frames = Parallel(n_jobs=n_jobs)(
                delayed(self._process_frame)(i, coords, elements, names)
                for i, coords in enumerate(tqdm(all_coords, 
                                               desc="Processing frames",
                                               unit="frame"))
            )
            
            # Write processed frames
            print("Writing trajectory...")
            with mda.Writer(str(output_path), n_atoms) as writer:
                # Create temporary universe for writing
                temp_u = mda.Universe.empty(n_atoms, 
                                           n_residues=1,
                                           atom_resindex=[0]*n_atoms,
                                           trajectory=True)
                temp_u.add_TopologyAttr('type', processed_frames[0][1])
                temp_u.add_TopologyAttr('name', processed_frames[0][2])
                
                for coords, elements, names in tqdm(processed_frames, 
                                                   desc="Writing frames", 
                                                   unit="frame"):
                    temp_u.atoms.positions = coords
                    temp_u.atoms.types = elements
                    temp_u.atoms.names = names
                    writer.write(temp_u.atoms)
        
        print(f"Trajectory written to {output_path}")


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
