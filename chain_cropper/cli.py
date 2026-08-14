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
from typing import List, Optional, Tuple

from .chain_cropper import TrajectoryProcessor
from .instrumentation import add_instrumentation_args, apply_instrumentation_args

try:
    import MDAnalysis as mda
except ImportError:
    print("Error: MDAnalysis is required but not installed.")
    print("Please install it with: pip install MDAnalysis")
    sys.exit(1)


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


def generate_output_paths(topology_path: Path, trajectory_path: Optional[Path],
                         output_arg: Optional[str]) -> Tuple[Path, Optional[Path]]:
    """
    Generate output file paths with _cropped suffix.
    
    Parameters
    ----------
    topology_path : Path
        Input topology file path
    trajectory_path : Optional[Path]
        Input trajectory file path
    output_arg : Optional[str]
        User-provided output path/directory
        
    Returns
    -------
    Tuple[Path, Optional[Path]]
        Output paths for (topology, trajectory)
    """
    if output_arg:
        output = Path(output_arg)
        
        # If output is a directory, generate filenames
        if output.is_dir() or (not output.exists() and not output.suffix):
            output.mkdir(parents=True, exist_ok=True)
            topo_out = output / f"{topology_path.stem}_cropped{topology_path.suffix}"
            traj_out = output / f"{trajectory_path.stem}_cropped{trajectory_path.suffix}" if trajectory_path else None
        else:
            # Output is a single file.
            if trajectory_path:
                # Both a topology and a trajectory output are needed from
                # one -o value -- always derive both from the INPUT
                # files' own extensions, not whatever extension -o
                # happens to have. Using output's extension verbatim for
                # the topology (e.g. `-o cropped.trr`) would silently
                # produce a trajectory-shaped "topology" file and no
                # actual cropped structure at all, leaving the cropped
                # trajectory with no matching topology to open it with.
                topo_out = output.parent / f"{output.stem}{topology_path.suffix}"
                traj_out = output.parent / f"{output.stem}_traj{trajectory_path.suffix}"
            else:
                # Single structure only -- honor the user's exact filename.
                topo_out = output
                traj_out = None
    else:
        # No output specified - use input directory with _cropped suffix
        topo_out = topology_path.parent / f"{topology_path.stem}_cropped{topology_path.suffix}"
        traj_out = trajectory_path.parent / f"{trajectory_path.stem}_cropped{trajectory_path.suffix}" if trajectory_path else None
    
    return topo_out, traj_out


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
  chain-cropper input.gro --chain-type alkyl --max-length 1
  
  # Process structure + trajectory (outputs both cropped files)
  chain-cropper input.gro input.xtc --chain-type ether --max-length 2
  
  # Specify output directory
  chain-cropper input.gro input.xtc -o results/
  
  # Remove all side chains (max-length 0)
  chain-cropper input.pdb --max-length 0
  
  # Batch process all files in current directory
  chain-cropper --batch -o results/
  
  # Verbose output
  chain-cropper input.gro --verbose
        """
    )
    
    # Input files
    parser.add_argument('topology', nargs='?', 
                       help='Topology file (gro, pdb, mol2, xyz)')
    parser.add_argument('trajectory', nargs='?', 
                       help='Trajectory file (xtc, trr, dcd) - optional')
    
    # Output
    parser.add_argument('-o', '--output', 
                       help='Output file or directory (default: same dir as input with _cropped suffix)')
    
    # Processing options
    parser.add_argument('--chain-type', choices=['alkyl', 'ether'], default='alkyl',
                       help='Type of chains to crop (default: alkyl)')
    parser.add_argument('--max-length', type=int, default=1,
                       help='Maximum number of heavy atoms to keep in chains (default: 1)')
    parser.add_argument('--cap-distance', type=float, default=1.09,
                       help='Distance for capping hydrogens in Angstrom (default: 1.09)')
    parser.add_argument('-j', '--n-jobs', type=int, default=-1,
                       help='Number of parallel jobs for trajectory processing. '
                            '-1 uses all cores, 1 disables parallelization (default: -1)')
    
    # Mode options
    parser.add_argument('--batch', action='store_true',
                       help='Process all topology files in current directory')

    add_instrumentation_args(parser)

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
                       output_arg: Optional[str], chain_type: str, max_length: int,
                       cap_distance: float, n_jobs: int) -> bool:
    """
    Process a single file or file pair.
    
    Parameters
    ----------
    topology_path : Path
        Path to topology file
    trajectory_path : Optional[Path]
        Path to trajectory file (optional)
    output_arg : Optional[str]
        Output path or directory
    chain_type : str
        Type of chains to crop
    max_length : int
        Maximum chain length
    cap_distance : float
        Capping distance
    n_jobs : int
        Number of parallel jobs
        
    Returns
    -------
    bool
        True if successful, False otherwise
    """
    logger = logging.getLogger(__name__)
    
    try:
        # Generate output paths
        topo_output, traj_output = generate_output_paths(
            topology_path, trajectory_path, output_arg
        )

        # Initialize processor
        processor = TrajectoryProcessor(cap_distance=cap_distance)

        # Analyse once, then write the structure. `process_trajectory` below
        # re-derives the same (structure, chain_type, max_length) key and
        # hits the cache, so the expensive analysis step never runs twice.
        logger.info(f"Processing structure with chain_type={chain_type}, max_length={max_length}")
        structure_universe = processor.analyze(
            structure_file=str(topology_path),
            chain_type=chain_type,
            max_chain_length=max_length,
        )
        processor.write_structure(structure_universe, topo_output)
        logger.info(f"Cropped structure written to {topo_output}")

        # Process trajectory if provided
        if trajectory_path:
            logger.info(f"Processing trajectory from {trajectory_path}")
            if n_jobs != 1:
                logger.info(f"Using {n_jobs if n_jobs > 0 else 'all available'} CPU cores for parallel processing")
            processor.process_trajectory(
                structure_file=str(topology_path),
                output_path=str(traj_output),
                trajectory_file=str(trajectory_path),
                chain_type=chain_type,
                max_chain_length=max_length,
                n_jobs=n_jobs,
                structure_universe=structure_universe,
            )
            logger.info(f"Cropped trajectory written to {traj_output}")

        processor.timing.log_summary()

        return True

    except Exception as e:
        logger.error(f"Error processing {topology_path}: {e}")
        if logger.level == logging.DEBUG:
            import traceback
            traceback.print_exc()
        return False


def process_batch(output_directory: Optional[str], chain_type: str, max_length: int,
                 cap_distance: float, n_jobs: int) -> int:
    """
    Process all topology files in current directory.
    
    Parameters
    ----------
    output_directory : Optional[str]
        Output directory (None uses current directory)
    chain_type : str
        Type of chains to crop
    max_length : int
        Maximum chain length
    cap_distance : float
        Capping distance
    n_jobs : int
        Number of parallel jobs
        
    Returns
    -------
    int
        Number of files processed successfully
    """
    logger = logging.getLogger(__name__)
    
    # Find topology files
    topology_files = get_topology_files(Path('.'))
    
    if not topology_files:
        logger.warning("No topology files found in current directory")
        return 0
    
    logger.info(f"Found {len(topology_files)} topology files")
    
    successful = 0
    for topo_file in topology_files:
        logger.info(f"Processing {topo_file.name}")
        
        if process_single_file(topo_file, None, output_directory, 
                              chain_type, max_length, cap_distance, n_jobs):
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
        # Batch mode - output should be a directory or None
        if args.output:
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
    
    # Validate n_jobs
    if args.n_jobs < -1 or args.n_jobs == 0:
        print("Error: n-jobs must be -1 (all cores) or a positive integer")
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
    
    # Set up logging -- file only, deliberately no console handler; the
    # console is reserved for the tqdm bars in TrajectoryProcessor.
    apply_instrumentation_args(args)
    logger = logging.getLogger(__name__)
    
    # Validate arguments
    try:
        validate_arguments(args)
    except SystemExit:
        return 1
    
    logger.info("Starting ChainCropper")
    
    try:
        if args.batch:
            # Batch processing mode
            successful = process_batch(args.output, args.chain_type, 
                                     args.max_length, args.cap_distance, args.n_jobs)
            
            if successful == 0:
                logger.error("No files processed successfully")
                return 1
            
        else:
            # Single file mode
            topology_path = Path(args.topology)
            trajectory_path = Path(args.trajectory) if args.trajectory else None
            
            success = process_single_file(topology_path, trajectory_path,
                                        args.output, args.chain_type,
                                        args.max_length, args.cap_distance, args.n_jobs)
            
            if not success:
                return 1
        
        logger.info("ChainCropper completed successfully")
        return 0
        
    except KeyboardInterrupt:
        logger.info("Process interrupted by user")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
