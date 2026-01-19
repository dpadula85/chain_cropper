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
import multiprocessing as mp
from typing import List, Optional, Tuple, Union, Dict, Set
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False


def _process_frame_worker(frame_idx: int, positions: np.ndarray, 
                          keep_indices: List[int], replace_info: List[Tuple[int, int]],
                          cap_distance: float, types: np.ndarray, 
                          names: np.ndarray) -> Tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    """
    Worker function to process a single frame.
    
    Parameters
    ----------
    frame_idx : int
        Frame index for ordering
    positions : np.ndarray
        Atomic positions for this frame
    keep_indices : List[int]
        Indices of atoms to keep
    replace_info : List[Tuple[int, int]]
        List of (anchor_idx, deleted_idx) pairs for capping
    cap_distance : float
        Distance for hydrogen caps
    types : np.ndarray
        Atom types array
    names : np.ndarray
        Atom names array
        
    Returns
    -------
    Tuple[int, np.ndarray, np.ndarray, np.ndarray]
        Frame index, positions, types, and names for kept atoms
    """
    coords = positions.copy()
    elements = types.copy()
    atom_names = names.copy()
    
    # Apply hydrogen capping
    for anchor_idx, deleted_idx in replace_info:
        anchor_pos = coords[anchor_idx]
        old_pos = coords[deleted_idx]
        
        # Calculate new hydrogen position
        bond_vector = old_pos - anchor_pos
        bond_vector = bond_vector / np.linalg.norm(bond_vector)
        new_h_pos = anchor_pos + bond_vector * cap_distance
        
        coords[deleted_idx] = new_h_pos
        elements[deleted_idx] = 'H'
        atom_names[deleted_idx] = 'H'
    
    # Extract only kept atoms
    final_coords = coords[keep_indices]
    final_types = elements[keep_indices]
    final_names = atom_names[keep_indices]
    
    return frame_idx, final_coords, final_types, final_names


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
            universe, chain_type, max_chain_length)
        
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
        
        return new_universe


class TrajectoryProcessor:
    """
    Handler for processing trajectories with chain cropping.
    Optimized with connectivity caching and parallel frame processing.
    """
    
    def __init__(self, cropper):
        """
        Initialize with a ChainCropper instance.
        
        Parameters
        ----------
        cropper : ChainCropper
            ChainCropper instance to use for processing
        """
        self.cropper = cropper
        self._cached_topology = None
    
    def process_trajectory(self, universe: mda.Universe, output_path: str,
                          chain_type: str = 'alkyl', max_chain_length: int = 1,
                          n_workers: int = None, show_progress: bool = True) -> None:
        """
        Process entire trajectory and write to file with parallelization.
        
        Parameters
        ----------
        universe : mda.Universe
            Input universe with trajectory
        output_path : str
            Output file path
        chain_type : str, default='alkyl'
            Type of chains to crop
        max_chain_length : int, default=1
            Maximum chain length to keep
        n_workers : int, optional
            Number of parallel workers. If None, uses CPU count - 1
        show_progress : bool, default=True
            Show progress bar during processing
        """
        # Determine output format
        output_path = Path(output_path)
        extension = output_path.suffix.lower()
        
        # Process first frame to determine topology (connectivity computed once)
        universe.trajectory[0]
        keep_indices, delete_indices, replace_indices = self.cropper.identify_chains_to_crop(
            universe, chain_type, max_chain_length)
        
        # Cache topology information
        self._cache_topology(universe, keep_indices, replace_indices, delete_indices)
        
        if extension in ['.xtc', '.trr', '.dcd']:
            self._write_trajectory_parallel(universe, output_path, n_workers, show_progress)
        elif extension in ['.xyz', '.gro', '.pdb']:
            self._write_single_frame(universe, output_path)
        else:
            raise ValueError(f"Unsupported output format: {extension}")
    
    def _cache_topology(self, universe: mda.Universe, keep_indices: List[int],
                       replace_indices: List[int], delete_indices: List[int]) -> None:
        """
        Cache topology information from first frame for reuse.
        
        Parameters
        ----------
        universe : mda.Universe
            Universe object
        keep_indices : List[int]
            Atoms to keep
        replace_indices : List[int]
            Atoms that need hydrogen caps
        delete_indices : List[int]
            Atoms to delete
        """
        # Build replacement information: (anchor_atom, deleted_atom) pairs
        replace_info = []
        atoms_to_keep = set(keep_indices)
        
        for atom_idx in replace_indices:
            connected = self.cropper.connectivity[atom_idx]
            connected_atoms = [idx for idx in connected if idx >= 0]
            
            deleted_connected = [idx for idx in connected_atoms if idx in delete_indices]
            if deleted_connected:
                replace_info.append((atom_idx, deleted_connected[0]))
        
        # Update keep_indices to include capped hydrogens
        final_keep_set = set(keep_indices)
        for _, deleted_idx in replace_info:
            final_keep_set.add(deleted_idx)
        
        final_keep_indices = sorted(list(final_keep_set))
        
        self._cached_topology = {
            'keep_indices': final_keep_indices,
            'replace_info': replace_info,
            'types': universe.atoms.types.copy(),
            'names': universe.atoms.names.copy()
        }
    
    def _write_trajectory_parallel(self, universe: mda.Universe, 
                                   output_path: Path, n_workers: int = None,
                                   show_progress: bool = True) -> None:
        """
        Write trajectory file using parallel processing.
        
        Parameters
        ----------
        universe : mda.Universe
            Input universe
        output_path : Path
            Output file path
        n_workers : int, optional
            Number of workers
        show_progress : bool, default=True
            Show progress bar
        """
        if n_workers is None:
            n_workers = max(1, mp.cpu_count() - 1)
        
        n_frames = len(universe.trajectory)
        keep_indices = self._cached_topology['keep_indices']
        replace_info = self._cached_topology['replace_info']
        cap_distance = self.cropper.cap_distance
        types = self._cached_topology['types']
        names = self._cached_topology['names']
        
        # Collect all positions from trajectory
        all_positions = []
        
        if show_progress and TQDM_AVAILABLE:
            pbar_load = tqdm(universe.trajectory, desc="Loading frames", unit="frame")
            for ts in pbar_load:
                all_positions.append(universe.atoms.positions.copy())
        else:
            for ts in universe.trajectory:
                all_positions.append(universe.atoms.positions.copy())
        
        # Process frames in parallel
        frame_results = [None] * n_frames
        
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {}
            for frame_idx in range(n_frames):
                future = executor.submit(
                    _process_frame_worker,
                    frame_idx,
                    all_positions[frame_idx],
                    keep_indices,
                    replace_info,
                    cap_distance,
                    types,
                    names
                )
                futures[future] = frame_idx
            
            # Collect results maintaining order with progress bar
            if show_progress and TQDM_AVAILABLE:
                pbar_process = tqdm(total=n_frames, desc="Processing frames", unit="frame")
                for future in as_completed(futures):
                    frame_idx, coords, frame_types, frame_names = future.result()
                    frame_results[frame_idx] = (coords, frame_types, frame_names)
                    pbar_process.update(1)
                pbar_process.close()
            else:
                for future in as_completed(futures):
                    frame_idx, coords, frame_types, frame_names = future.result()
                    frame_results[frame_idx] = (coords, frame_types, frame_names)
        
        # Write trajectory in order
        n_atoms = len(keep_indices)
        
        with mda.Writer(str(output_path), n_atoms) as writer:
            # Create temporary universe for writing
            temp_u = mda.Universe.empty(n_atoms, trajectory=True)
            temp_u.add_TopologyAttr('type', frame_results[0][1])
            temp_u.add_TopologyAttr('name', frame_results[0][2])
            
            if show_progress and TQDM_AVAILABLE:
                pbar_write = tqdm(frame_results, desc="Writing trajectory", unit="frame")
                for coords, frame_types, frame_names in pbar_write:
                    temp_u.atoms.positions = coords
                    writer.write(temp_u.atoms)
            else:
                for coords, frame_types, frame_names in frame_results:
                    temp_u.atoms.positions = coords
                    writer.write(temp_u.atoms)
    
    def _write_single_frame(self, universe: mda.Universe, output_path: Path) -> None:
        """Write single frame file using cached topology."""
        keep_indices = self._cached_topology['keep_indices']
        replace_info = self._cached_topology['replace_info']
        
        coords = universe.atoms.positions.copy()
        elements = self._cached_topology['types'].copy()
        atom_names = self._cached_topology['names'].copy()
        
        # Apply hydrogen capping
        for anchor_idx, deleted_idx in replace_info:
            anchor_pos = coords[anchor_idx]
            old_pos = coords[deleted_idx]
            
            bond_vector = old_pos - anchor_pos
            bond_vector = bond_vector / np.linalg.norm(bond_vector)
            new_h_pos = anchor_pos + bond_vector * self.cropper.cap_distance
            
            coords[deleted_idx] = new_h_pos
            elements[deleted_idx] = 'H'
            atom_names[deleted_idx] = 'H'
        
        # Create new universe with kept atoms
        temp_u = mda.Universe.empty(len(keep_indices), trajectory=True)
        temp_u.add_TopologyAttr('type', elements[keep_indices])
        temp_u.add_TopologyAttr('name', atom_names[keep_indices])
        temp_u.atoms.positions = coords[keep_indices]
        
        temp_u.atoms.write(str(output_path))
    
if __name__ == '__main__':
    pass
