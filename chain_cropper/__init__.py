"""
ChainCropper - A tool for cropping molecular side chains from structures and trajectories.

This package provides functionality to crop alkyl and ether side chains from molecular
structures while properly capping broken bonds with hydrogen atoms. It supports both
single structures and molecular dynamics trajectories.

Classes:
    ChainCropper: Main class for cropping side chains
    TrajectoryProcessor: Handler for processing trajectories
"""

__version__ = "0.1"
__author__ = "Daniele Padula"
__email__ = "daniele.padula@unisi.it"

from .chain_cropper import ChainCropper, TrajectoryProcessor

__all__ = ["ChainCropper", "TrajectoryProcessor"]
