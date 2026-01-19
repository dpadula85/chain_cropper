#!/usr/bin/env python3
"""
Command Line Interface for ChainCropper.

This module provides a command-line interface for cropping molecular side chains
from structures and trajectories.
"""

import sys
import logging
import argparse
from pathlib import Path
from typing import List, Optional

from .chain_cropper import ChainCropper, TrajectoryProcessor

try:
    import MDAnalysis as mda
except ImportError:
    print("Error: MDAnalysis is required but not installed.")
    print("Please install it with: pip install MDAnalysis")
    sys.exit(1)


def setup_logging(verbose: bool = False) -> None:
    """
    Set up logging configuration.
    
    Parameters
    ----------
    verbose : bool
        Enable verbose logging
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def validate_input_file(filepath: str) -> Path:
    """
    Validate that input file exists and has supported extension.
    
    Parameters
    ----------
    filepath : str
        Path to input file
        
    Returns
    -------
    Path
        Validated Path object
        
    Raises
    ------
    FileNotFoundError
        If file doesn't exist
    ValueError
        If file extension is not supported
    """
    path = Path(filepath)
    
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {filepath}")
    
    supported_topo_exts = {'.gro', '.pdb', '.mol2', '.xyz'}
    supported_traj_exts = {'.xtc', '.trr', '.dcd'}
    
    if path.suffix.lower() not in supported_topo_exts.union(supported_traj_exts):
        raise ValueError(f"Unsupported file format: {path.suffix}")
    
    return path


def get_topology_files(directory: Path) -> List[Path]:
    """
    Get all topology files in a directory.
    
    Parameters
    ----------
    directory : Path
        Directory to search
        
    Returns
    -------
    List[Path]
        List of topology files
    """
    topology_extensions = ['.gro', '.pdb', '.mol2', '.xyz']
    topology_files = []
    
    for ext in topology_extensions:
        topology_files.extend(directory.glob(f'*{ext}'))
    
    return sorted(topology_files)


def create_cli_parser() -> argparse.ArgumentParser:
    """
    Create command line interface parser.
    
    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser
    """
    parser = argparse.ArgumentParser(
        description="Crop molecular side chains from structures and trajectories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Crop alkyl chains, keeping 1 carbon
  chain-cropper input.gro -o output.gro --chain-type alkyl --max-length 1
  
  # Process trajectory with ether chain cropping
  chain-cropper input.gro input.xtc -o output.xtc --chain-type ether --max-length 2
  
  # Remove all side chains (max-length 0)
  chain-cropper input.pdb -o output.xyz --max-length 0
  
  # Process trajectory with 8 parallel workers
  chain-cropper input.gro input.xtc -o output.xtc --workers 8
  
  # Disable progress bar
  chain-cropper input.gro input.xtc -o output.xtc --no-progress
  
  # Batch process all files in current directory
  chain-cropper --batch -o results/
  
  # Verbose output
  chain-cropper input.gro -o output.gro --verbose
        """
    )
    
    # Input files
    parser.add_argument('topology', nargs='?', 
                       help='Topology file (gro, pdb, mol2, xyz)')
    parser.add_argument('trajectory', nargs='?', 
                       help='Trajectory file (xtc, trr, dcd) - optional')
    
    # Output
    parser.add_argument('-o', '--output', required=True,
                       help='Output file path or directory (for batch mode)')
    
    # Processing options
    parser.add_argument('--chain-type', choices=['alkyl', 'ether'], default='alkyl',
                       help='Type of chains to crop (default: alkyl)')
    parser.add_argument('--max-length', type=int, default=1,
                       help='Maximum number of heavy atoms to keep in chains (default: 1)')
    parser.add_argument('--cap-distance', type=float, default=1.09,
                       help='Distance for capping hydrogens in Angstrom (default: 1.09)')
    
    # Performance options
    parser.add_argument('--workers', type=int, default=None,
                       help='Number of parallel workers for trajectory processing (default: CPU count - 1)')
    parser.add_argument('--no-progress', action='store_true',
                       help='Disable progress bar')
    
    # Mode options
    parser.add_argument('--batch', action='store_true',
                       help='Process all topology files in current directory')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose output')
    
    # Version
    try:
        from . import __version__
        parser.add_argument('--version', action='version', 
                           version=f'ChainCropper {__version__}')
    except ImportError:
        parser.add_argument('--version', action='version', 
                           version='ChainCropper 1.0.0')
    
    return parser


def process_single_file(topology_path: Path, trajectory_path: Optional[Path],
                       output_path: Path, chain_type: str, max_length: int,
                       cap_distance: float, n_workers: Optional[int],
                       show_progress: bool) -> bool:
    """
    Process a single file or file pair.
    
    Parameters
    ----------
    topology_path : Path
        Path to topology file
    trajectory_path : Optional[Path]
        Path to trajectory file (optional)
    output_path : Path
        Output file path
    chain_type : str
        Type of chains to crop
    max_length : int
        Maximum chain length
    cap_distance : float
        Capping distance
    n_workers : Optional[int]
        Number of parallel workers
    show_progress : bool
        Show progress bar
        
    Returns
    -------
    bool
        True if successful, False otherwise
    """
    logger = logging.getLogger(__name__)
    
    try:
        logger.info(f"Loading universe from {topology_path}")
        
        # Load universe
        if trajectory_path:
            logger.info(f"Including trajectory: {trajectory_path}")
            universe = mda.Universe(str(topology_path), str(trajectory_path))
        else:
            universe = mda.Universe(str(topology_path))
        
        logger.info(f"Universe loaded: {len(universe.atoms)} atoms, "
                   f"{len(universe.trajectory)} frames")
        
        # Initialize processor
        cropper = ChainCropper(cap_distance=cap_distance)
        processor = TrajectoryProcessor(cropper)
        
        # Process
        if n_workers is not None:
            logger.info(f"Using {n_workers} parallel workers")
        
        logger.info(f"Processing with chain_type={chain_type}, "
                   f"max_length={max_length}")
        
        processor.process_trajectory(universe, str(output_path), 
                                   chain_type, max_length,
                                   n_workers=n_workers,
                                   show_progress=show_progress)
        
        logger.info(f"Successfully written to {output_path}")
        return True
        
    except Exception as e:
        logger.error(f"Error processing {topology_path}: {e}")
        return False


def process_batch(output_directory: Path, chain_type: str, max_length: int,
                 cap_distance: float, n_workers: Optional[int],
                 show_progress: bool) -> int:
    """
    Process all topology files in current directory.
    
    Parameters
    ----------
    output_directory : Path
        Output directory
    chain_type : str
        Type of chains to crop
    max_length : int
        Maximum chain length
    cap_distance : float
        Capping distance
    n_workers : Optional[int]
        Number of parallel workers
    show_progress : bool
        Show progress bar
        
    Returns
    -------
    int
        Number of files processed successfully
    """
    logger = logging.getLogger(__name__)
    
    # Create output directory if it doesn't exist
    output_directory.mkdir(parents=True, exist_ok=True)
    
    # Find topology files
    topology_files = get_topology_files(Path('.'))
    
    if not topology_files:
        logger.warning("No topology files found in current directory")
        return 0
    
    logger.info(f"Found {len(topology_files)} topology files")
    
    successful = 0
    for topo_file in topology_files:
        # Generate output filename
        output_file = f"cropped_{topo_file.stem}.gro"
        output_path = output_directory / output_file
        
        logger.info(f"Processing {topo_file.name} -> {output_file}")
        
        if process_single_file(topo_file, None, output_path, 
                              chain_type, max_length, cap_distance,
                              n_workers, show_progress):
            successful += 1
    
    logger.info(f"Batch processing complete: {successful}/{len(topology_files)} successful")
    return successful


def validate_arguments(args: argparse.Namespace) -> None:
    """
    Validate command line arguments.
    
    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments
        
    Raises
    ------
    SystemExit
        If arguments are invalid
    """
    if args.batch:
        # Batch mode - output should be a directory
        output_path = Path(args.output)
        if output_path.exists() and not output_path.is_dir():
            print("Error: In batch mode, output must be a directory")
            sys.exit(1)
    else:
        # Single file mode - need topology
        if not args.topology:
            print("Error: Topology file required when not in batch mode")
            sys.exit(1)
        
        # Validate input files
        try:
            validate_input_file(args.topology)
            if args.trajectory:
                validate_input_file(args.trajectory)
        except (FileNotFoundError, ValueError) as e:
            print(f"Error: {e}")
            sys.exit(1)
    
    # Validate max_length
    if args.max_length < 0:
        print("Error: max-length must be non-negative")
        sys.exit(1)
    
    # Validate cap_distance
    if args.cap_distance <= 0:
        print("Error: cap-distance must be positive")
        sys.exit(1)
    
    # Validate workers
    if args.workers is not None and args.workers < 1:
        print("Error: workers must be at least 1")
        sys.exit(1)


def main() -> int:
    """
    Main function for CLI usage.
    
    Returns
    -------
    int
        Exit code (0 for success, 1 for error)
    """
    parser = create_cli_parser()
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)
    
    # Validate arguments
    try:
        validate_arguments(args)
    except SystemExit:
        return 1
    
    logger.info("Starting ChainCropper")
    
    # Determine progress bar setting
    show_progress = not args.no_progress
    
    try:
        if args.batch:
            # Batch processing mode
            output_directory = Path(args.output)
            successful = process_batch(output_directory, args.chain_type, 
                                     args.max_length, args.cap_distance,
                                     args.workers, show_progress)
            
            if successful == 0:
                logger.error("No files processed successfully")
                return 1
            
        else:
            # Single file mode
            topology_path = Path(args.topology)
            trajectory_path = Path(args.trajectory) if args.trajectory else None
            output_path = Path(args.output)
            
            success = process_single_file(topology_path, trajectory_path,
                                        output_path, args.chain_type,
                                        args.max_length, args.cap_distance,
                                        args.workers, show_progress)
            
            if not success:
                return 1
        
        logger.info("ChainCropper completed successfully")
        return 0
        
    except KeyboardInterrupt:
        logger.info("Process interrupted by user")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
