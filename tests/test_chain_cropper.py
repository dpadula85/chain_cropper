"""
Regression tests for ChainCropper/TrajectoryProcessor's core cropping
logic, on a small hand-verified synthetic molecule:

    N(0) - C1(1) - C2(2) - C3(5), each chain carbon also bonded to H's
           (3,4)          (6,7)   (8,9,10)

With chain_type='alkyl', max_chain_length=1: C1 is kept (it's the first
atom of the alkyl chain off the N anchor), C2/C3 and their hydrogens are
deleted, and C1 gets a new capping hydrogen at `cap_distance` where C2
used to be (repurposing C2's own atom index, not appending a new atom --
see `TrajectoryProcessor.final_indices`'s docstring).

Written for the identify_chains_to_crop/_apply_cropping performance fix
(hoisting `universe.atoms.types` out of the atoms_to_delete loop, using a
set instead of a list for the delete-index membership check, and direct
array indexing instead of a giant selection string) -- these assert the
fix did not change behaviour, not just that it's faster.
"""
import numpy as np
import pytest
import MDAnalysis as mda

from chain_cropper.chain_cropper import TrajectoryProcessor

CAP_DISTANCE = 1.09


@pytest.fixture
def universe():
    n_atoms = 11
    u = mda.Universe.empty(n_atoms, trajectory=True)
    types = ["N", "C", "C", "H", "H", "C", "H", "H", "H", "H", "H"]
    u.add_TopologyAttr("type", types)
    u.add_TopologyAttr("name", types)
    u.atoms.positions = np.random.RandomState(0).rand(n_atoms, 3) * 3
    u.add_bonds([(0, 1), (1, 2), (1, 3), (1, 4), (2, 5), (2, 6), (2, 7),
                 (5, 8), (5, 9), (5, 10)])
    return u


def test_identify_chains_to_crop(universe):
    processor = TrajectoryProcessor(cap_distance=CAP_DISTANCE)
    keep, delete, replace = processor.identify_chains_to_crop(universe, "alkyl", 1)

    assert sorted(int(i) for i in keep) == [0, 1, 3, 4]
    assert sorted(int(i) for i in delete) == [2, 5, 6, 7, 8, 9, 10]
    assert sorted(int(i) for i in replace) == [1]


def test_apply_cropping_caps_at_configured_distance(universe):
    processor = TrajectoryProcessor(cap_distance=CAP_DISTANCE)
    keep, delete, replace = processor.identify_chains_to_crop(universe, "alkyl", 1)
    processor.keep_indices = keep
    processor.delete_indices = delete
    processor.replace_indices = replace

    original_c1_pos = universe.atoms.positions[1].copy()
    cropped = processor._apply_cropping(universe)

    assert cropped.atoms.n_atoms == 5
    assert list(cropped.atoms.names) == ["N", "C", "H", "H", "H"]
    assert list(cropped.atoms.elements) == ["N", "C", "H", "H", "H"]

    # The capped H repurposes original atom 2 (C2)'s index -- final_indices
    # records that, so its position in the cropped structure is found by
    # looking up where original index 2 landed.
    cap_pos_in_cropped = list(processor.final_indices).index(2)
    cap_pos = cropped.atoms.positions[cap_pos_in_cropped]
    assert np.linalg.norm(cap_pos - original_c1_pos) == pytest.approx(CAP_DISTANCE, rel=1e-5)


def test_apply_cropping_preserves_box(universe):
    universe.dimensions = [30.0, 30.0, 30.0, 90.0, 90.0, 90.0]
    processor = TrajectoryProcessor(cap_distance=CAP_DISTANCE)
    keep, delete, replace = processor.identify_chains_to_crop(universe, "alkyl", 1)
    processor.keep_indices, processor.delete_indices, processor.replace_indices = keep, delete, replace

    cropped = processor._apply_cropping(universe)
    assert cropped.dimensions[:3] == pytest.approx([30.0, 30.0, 30.0])
