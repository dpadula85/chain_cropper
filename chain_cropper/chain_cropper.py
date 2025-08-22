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
    """
    
    def __init__(self, cropper: ChainCropper):
        """
        Initialize with a ChainCropper instance.
        
        Parameters
        ----------
        cropper : ChainCropper
            ChainCropper instance to use for processing
        """
        self.cropper = cropper
    
    def process_trajectory(self, universe: mda.Universe, output_path: str,
                          chain_type: str = 'alkyl', max_chain_length: int = 1) -> None:
        """
        Process entire trajectory and write to file.
        
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
        """
        # Process first frame to determine structure
        universe.trajectory[0]
        cropped_universe = self.cropper.crop_chains(universe, chain_type, max_chain_length)
        
        # Determine output format
        output_path = Path(output_path)
        extension = output_path.suffix.lower()
        
        if extension in ['.xtc', '.trr', '.dcd']:
            self._write_trajectory(universe, cropped_universe, output_path, 
                                 chain_type, max_chain_length)
        elif extension in ['.xyz', '.gro', '.pdb']:
            self._write_single_frame(cropped_universe, output_path)
        else:
            raise ValueError(f"Unsupported output format: {extension}")
    
    def _write_trajectory(self, original_universe: mda.Universe, 
                         template_universe: mda.Universe, output_path: Path,
                         chain_type: str, max_chain_length: int) -> None:
        """Write trajectory file."""
        n_atoms = len(template_universe.atoms)
        
        with mda.Writer(str(output_path), n_atoms) as writer:
            for ts in original_universe.trajectory:
                # Process current frame
                frame_universe = self.cropper.crop_chains(
                    original_universe, chain_type, max_chain_length)
                writer.write(frame_universe.atoms)
    
    def _write_single_frame(self, universe: mda.Universe, output_path: Path) -> None:
        """Write single frame file."""
        extension = output_path.suffix.lower()
        universe.atoms.write(str(output_path))
    
if __name__ == '__main__':
    pass
