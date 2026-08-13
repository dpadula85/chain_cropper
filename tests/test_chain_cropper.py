"""
Regression tests for ChainCropper/TrajectoryProcessor's core cropping logic.

Most of these pin down behaviour that was previously wrong. The molecules
are small and hand-verified; each test's docstring says which failure it
guards against.
"""
import numpy as np
import pytest
import MDAnalysis as mda

from chain_cropper.chain_cropper import ChainCropper, TrajectoryProcessor
from chain_cropper.topology import build_connectivity, get_sp2, side_chain_atoms

CAP_DISTANCE = 1.09


def make_universe(types, bonds, seed=0, box=None):
    """Universe with explicit bonds and arbitrary but reproducible geometry."""
    n_atoms = len(types)
    u = mda.Universe.empty(n_atoms, trajectory=True)
    u.add_TopologyAttr("type", list(types))
    u.add_TopologyAttr("name", list(types))
    u.atoms.positions = np.random.RandomState(seed).rand(n_atoms, 3) * 10
    u.add_bonds(bonds)
    if box is not None:
        u.dimensions = box
    return u


# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------

@pytest.fixture
def propyl():
    """
    N(0) - C1(1) - C2(2) - C3(5), each chain carbon carrying hydrogens:

        N(0) - C1(1)(H3,H4) - C2(2)(H6,H7) - C3(5)(H8,H9,H10)

    A single unbranched alkyl chain hanging off a non-chain anchor.
    """
    types = ["N", "C", "C", "H", "H", "C", "H", "H", "H", "H", "H"]
    bonds = [(0, 1), (1, 2), (1, 3), (1, 4), (2, 5), (2, 6), (2, 7),
             (5, 8), (5, 9), (5, 10)]
    return make_universe(types, bonds)


def _aromatic_ring(first_c, first_h):
    """Six-membered all-carbon ring, five of its atoms carrying one H."""
    ring = [(first_c + i, first_c + (i + 1) % 6) for i in range(6)]
    hydrogens = [(first_c + i, first_h + i - 1) for i in range(1, 6)]
    return ring, hydrogens


def fluorene(bridge_substituents=0, arm_length=0):
    """
    Fluorene-like: two benzene rings fused through one four-coordinate
    bridge carbon (atom 12), the way 9,9-dialkylfluorene, CPDT and IDT are
    built. `bridge_substituents` alkyl arms of `arm_length` carbons each are
    hung off the bridge; the remaining bridge valences carry hydrogen.

    Returns (universe, bridge_index).
    """
    types = ["C"] * 13
    bonds = []
    for ring_start in (0, 6):
        ring, _ = _aromatic_ring(ring_start, 0)
        bonds += ring
    bonds += [(0, 12), (6, 12)]

    # One H on each aromatic carbon except the two fused to the bridge.
    for aromatic in [1, 2, 3, 4, 5, 7, 8, 9, 10, 11]:
        bonds.append((aromatic, len(types)))
        types.append("H")

    def add_hydrogens(carbon, count):
        for _ in range(count):
            bonds.append((carbon, len(types)))
            types.append("H")

    # Alkyl arms off the bridge: a chain of `arm_length` carbons, each
    # saturated with hydrogen (two for the internal ones, three for the tip).
    for _ in range(bridge_substituents):
        previous = 12
        arm = []
        for _ in range(arm_length):
            carbon = len(types)
            types.append("C")
            bonds.append((previous, carbon))
            arm.append(carbon)
            previous = carbon
        for position, carbon in enumerate(arm):
            add_hydrogens(carbon, 3 if position == len(arm) - 1 else 2)

    # Whatever bridge valences are left get hydrogen.
    add_hydrogens(12, 2 - bridge_substituents)

    return make_universe(types, bonds), 12


# --------------------------------------------------------------------------
# Connectivity
# --------------------------------------------------------------------------

def test_connectivity_keeps_more_than_four_bonds():
    """
    A five-coordinate atom used to lose its fifth bond without a word: the
    neighbour matrix was hard-coded four columns wide and the fill loop just
    `break`ed when it found no free slot.
    """
    u = make_universe(["S", "O", "O", "C", "C", "C"],
                      [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)])
    connectivity, degree = build_connectivity(u)

    assert degree[0] == 5
    assert sorted(int(i) for i in connectivity[0] if i >= 0) == [1, 2, 3, 4, 5]


def test_connectivity_ignores_duplicate_bonds():
    """A bond listed twice must not consume two neighbour slots."""
    u = make_universe(["C", "H"], [(0, 1), (0, 1), (1, 0)])
    connectivity, degree = build_connectivity(u)

    assert degree[0] == 1
    assert sorted(int(i) for i in connectivity[0] if i >= 0) == [1]


def test_connectivity_pads_to_at_least_four_columns():
    """Callers may still assume four columns when nothing needs more."""
    u = make_universe(["C", "H"], [(0, 1)])
    connectivity, _degree = build_connectivity(u)

    assert connectivity.shape == (2, 4)


# --------------------------------------------------------------------------
# sp2/sp3 classification and the ether branch
# --------------------------------------------------------------------------

def test_get_sp2_finds_saturated_and_unsaturated(propyl):
    """Chain carbons are saturated; the three-coordinate N is not."""
    sp2, sp3 = get_sp2(propyl)

    assert sorted(int(i) for i in sp3) == [1, 2, 5]
    assert sorted(int(i) for i in sp2) == [0]


def test_get_sp2_ether_branch_is_reachable():
    """
    `get_sp2(u, ether=True)` used to be a no-op: the copies of this function
    in SelIntCoords and oligomer_builder tested `if alkyl: ... elif ether:`,
    and `alkyl` defaults to True, so the ether branch could not be reached
    without also passing `alkyl=False`. No caller ever did.

    Anisole-like: aromatic C(0) - O(1) - C(2)H3. In ether mode the oxygen
    counts as side chain, so it must not be reported as sp2.
    """
    types = ["C", "O", "C", "H", "H", "H", "C", "C"]
    bonds = [(0, 1), (1, 2), (2, 3), (2, 4), (2, 5), (0, 6), (0, 7)]
    u = make_universe(types, bonds)

    sp2_alkyl, _ = get_sp2(u, ether=False)
    sp2_ether, _ = get_sp2(u, ether=True)

    assert 1 in set(int(i) for i in sp2_alkyl)
    assert 1 not in set(int(i) for i in sp2_ether)


def test_ether_mode_walks_through_oxygen():
    """
    Ar-O-CH2-CH3: in alkyl mode the walk cannot cross the oxygen, so the
    ethyl is two separate depth-1 stubs. In ether mode the oxygen is depth 1
    and the carbons are depth 2 and 3.
    """
    types = ["C", "C", "C", "O", "C", "C"] + ["H"] * 5
    bonds = [(0, 1), (1, 2), (2, 0),          # tiny unsaturated core
             (0, 3), (3, 4), (4, 5),          # -O-CH2-CH3
             (4, 6), (4, 7), (5, 8), (5, 9), (5, 10)]
    u = make_universe(types, bonds)

    alkyl = set(int(i) for i in side_chain_atoms(u, "alkyl"))
    ether = set(int(i) for i in side_chain_atoms(u, "ether"))

    assert 3 not in alkyl
    assert 3 in ether


def test_unknown_chain_type_raises():
    """It used to fall through to 'alkyl' for any unrecognised value."""
    u = make_universe(["C", "H"], [(0, 1)])
    with pytest.raises(ValueError, match="chain_type"):
        side_chain_atoms(u, "aromatic")


# --------------------------------------------------------------------------
# identify_chains_to_crop
# --------------------------------------------------------------------------

def test_identify_chains_to_crop(propyl):
    """C1 is kept as the one allowed chain carbon; C2/C3 and their H's go."""
    keep, delete, replace = TrajectoryProcessor().identify_chains_to_crop(
        propyl, "alkyl", 1)

    assert sorted(int(i) for i in keep) == [0, 1, 3, 4]
    assert sorted(int(i) for i in delete) == [2, 5, 6, 7, 8, 9, 10]
    assert sorted(int(i) for i in replace) == [1]


@pytest.mark.parametrize("max_chain_length", [0, 1, 2, 3, 4])
def test_every_atom_is_either_kept_or_deleted(propyl, max_chain_length):
    """
    The old walk classified an atom as neither when its chain was short
    enough to keep whole -- it was excluded from the keep set for being a
    chain atom, but nothing ever added it to the delete set. Such atoms
    disappeared from the output silently, hydrogens and all.
    """
    keep, delete, _replace = TrajectoryProcessor().identify_chains_to_crop(
        propyl, "alkyl", max_chain_length)

    assert set(keep) | set(delete) == set(range(propyl.atoms.n_atoms))
    assert not set(keep) & set(delete)


def test_max_chain_length_zero_caps_the_anchor(propyl):
    """
    `--max-length 0` is documented as "remove all side chains". It used to
    delete the chain and then cap nothing at all, because the capping loop
    only ever looked at kept *chain* atoms and at max_chain_length=0 there
    are none -- leaving the anchor with a dangling valence.
    """
    processor = TrajectoryProcessor(cap_distance=CAP_DISTANCE)
    keep, delete, replace = processor.identify_chains_to_crop(propyl, "alkyl", 0)

    assert sorted(int(i) for i in keep) == [0]
    assert sorted(int(i) for i in replace) == [0]

    anchor_position = propyl.atoms.positions[0].copy()
    processor.keep_indices, processor.delete_indices, processor.replace_indices = \
        keep, delete, replace
    cropped = processor._apply_cropping(propyl)

    assert cropped.atoms.n_atoms == 2
    assert sorted(cropped.atoms.names) == ["H", "N"]
    cap = cropped.atoms.positions[list(cropped.atoms.names).index("H")]
    assert np.linalg.norm(cap - anchor_position) == pytest.approx(
        CAP_DISTANCE, rel=1e-5)


def test_negative_max_chain_length_raises(propyl):
    with pytest.raises(ValueError, match="non-negative"):
        TrajectoryProcessor().identify_chains_to_crop(propyl, "alkyl", -1)


def test_branched_chain_deletes_both_branches():
    """
    Ar-CH2-CH(CH3)(CH3): the old walk followed exactly one branch out of
    each atom, so the branch it did not follow was left classified as
    neither kept nor deleted -- its carbon vanished from the output while
    its three hydrogens stayed behind, bonded to nothing.
    """
    types = ["C", "C", "C",                  # unsaturated core ring
             "C", "C", "C", "C"] + ["H"] * 9
    bonds = [(0, 1), (1, 2), (2, 0),
             (0, 3), (3, 4), (4, 5), (4, 6),
             (3, 7), (3, 8), (4, 9),
             (5, 10), (5, 11), (5, 12),
             (6, 13), (6, 14), (6, 15)]
    u = make_universe(types, bonds)

    keep, delete, _replace = TrajectoryProcessor().identify_chains_to_crop(
        u, "alkyl", 1)

    assert set(keep) | set(delete) == set(range(u.atoms.n_atoms))
    # Both terminal methyls, and every one of their hydrogens, are deleted.
    for atom in [4, 5, 6, 9, 10, 11, 12, 13, 14, 15]:
        assert atom in set(delete), f"atom {atom} should be deleted"


# --------------------------------------------------------------------------
# Bridge carbons (fluorene / CPDT / IDT)
# --------------------------------------------------------------------------

@pytest.mark.parametrize("max_chain_length", [0, 1, 2])
def test_fused_bridge_carbon_is_never_cropped(max_chain_length):
    """
    Fluorene's bridge carbon is four-coordinate, so a walk that treats every
    four-coordinate atom as side chain reaches it and either deletes it --
    splitting the conjugated unit into two uncapped halves -- or drops it
    silently. It holds two pieces of core together, so it is core.
    """
    u, bridge = fluorene(bridge_substituents=0)
    keep, delete, _replace = TrajectoryProcessor().identify_chains_to_crop(
        u, "alkyl", max_chain_length)

    assert bridge in set(keep)
    assert bridge not in set(delete)
    assert set(keep) | set(delete) == set(range(u.atoms.n_atoms))


def test_bridge_carbon_gets_one_cap_per_severed_bond():
    """
    9,9-diethylfluorene at max_chain_length=0: the bridge loses BOTH arms
    and so needs two capping hydrogens. The old code capped only the first
    severed bond per atom (`deleted_connected[0]`), leaving every
    gem-disubstituted side-chain carbon one hydrogen short.
    """
    u, bridge = fluorene(bridge_substituents=2, arm_length=2)
    processor = TrajectoryProcessor(cap_distance=CAP_DISTANCE)
    keep, delete, replace = processor.identify_chains_to_crop(u, "alkyl", 0)

    assert bridge in set(replace)
    bridge_position = u.atoms.positions[bridge].copy()

    processor.keep_indices, processor.delete_indices, processor.replace_indices = \
        keep, delete, replace
    cropped = processor._apply_cropping(u)

    # 12 aromatic C + 10 aromatic H + the bridge itself + 2 caps. The bridge
    # is disubstituted, so it carries no hydrogen of its own.
    assert cropped.atoms.n_atoms == 25

    cap_distances = sorted(
        np.linalg.norm(cropped.atoms.positions[i] - bridge_position)
        for i, name in enumerate(cropped.atoms.names) if name == "H"
    )
    # Both severed arms became a hydrogen at cap_distance. Capping only the
    # first would leave exactly one here, and the bridge one H short.
    assert cap_distances[:2] == pytest.approx([CAP_DISTANCE] * 2, rel=1e-5)
    assert cap_distances[2] > CAP_DISTANCE * 1.5


def test_cutting_a_saturated_ring_mid_ring_raises():
    """
    A cyclohexyl substituent cut so that one deleted carbon sits between two
    kept ones cannot be capped: each cap is carried by the deleted atom's own
    index, so one deleted atom cannot become two hydrogens. Better to say so
    than to hand back a structure quietly one hydrogen short.
    """
    types = ["C", "C", "C"] + ["C"] * 6 + ["H"] * 11
    ring = [(3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 3)]
    bonds = [(0, 1), (1, 2), (2, 0), (0, 3)] + ring
    n = 9
    for carbon, n_h in [(3, 1), (4, 2), (5, 2), (6, 2), (7, 2), (8, 2)]:
        for _ in range(n_h):
            bonds.append((carbon, n))
            n += 1
    u = make_universe(types, bonds)

    with pytest.raises(ValueError, match="between two kept atoms"):
        TrajectoryProcessor().identify_chains_to_crop(u, "alkyl", 3)


# --------------------------------------------------------------------------
# Building the cropped universe
# --------------------------------------------------------------------------

def test_apply_cropping_caps_at_configured_distance(propyl):
    processor = TrajectoryProcessor(cap_distance=CAP_DISTANCE)
    keep, delete, replace = processor.identify_chains_to_crop(propyl, "alkyl", 1)
    processor.keep_indices, processor.delete_indices, processor.replace_indices = \
        keep, delete, replace

    original_c1_position = propyl.atoms.positions[1].copy()
    cropped = processor._apply_cropping(propyl)

    assert cropped.atoms.n_atoms == 5
    assert list(cropped.atoms.names) == ["N", "C", "H", "H", "H"]
    assert list(cropped.atoms.elements) == ["N", "C", "H", "H", "H"]

    # The capping H repurposes original atom 2 (C2)'s index -- final_indices
    # records that, so its position is found by looking up where original
    # index 2 landed.
    cap = cropped.atoms.positions[list(processor.final_indices).index(2)]
    assert np.linalg.norm(cap - original_c1_position) == pytest.approx(
        CAP_DISTANCE, rel=1e-5)


def test_apply_cropping_preserves_box(propyl):
    propyl.dimensions = [30.0, 30.0, 30.0, 90.0, 90.0, 90.0]
    processor = TrajectoryProcessor(cap_distance=CAP_DISTANCE)
    keep, delete, replace = processor.identify_chains_to_crop(propyl, "alkyl", 1)
    processor.keep_indices, processor.delete_indices, processor.replace_indices = \
        keep, delete, replace

    cropped = processor._apply_cropping(propyl)
    assert cropped.dimensions[:3] == pytest.approx([30.0, 30.0, 30.0])


def test_final_indices_are_ascending_original_indices(propyl):
    processor = TrajectoryProcessor()
    keep, delete, replace = processor.identify_chains_to_crop(propyl, "alkyl", 1)
    processor.keep_indices, processor.delete_indices, processor.replace_indices = \
        keep, delete, replace
    cropped = processor._apply_cropping(propyl)

    indices = np.asarray(processor.final_indices)
    assert indices.tolist() == sorted(indices.tolist())
    assert len(indices) == cropped.atoms.n_atoms
    assert indices.min() >= 0
    assert indices.max() < propyl.atoms.n_atoms


def test_crop_chains_does_not_mutate_its_input(propyl):
    """
    It used to write the capped positions, types and names straight back
    onto the input universe, so a caller that cropped and then went back to
    the original found atom 2 turned into a hydrogen.
    """
    types_before = propyl.atoms.types.copy()
    names_before = propyl.atoms.names.copy()
    positions_before = propyl.atoms.positions.copy()

    ChainCropper(cap_distance=CAP_DISTANCE).crop_chains(propyl, "alkyl", 1)

    assert list(propyl.atoms.types) == list(types_before)
    assert list(propyl.atoms.names) == list(names_before)
    assert propyl.atoms.positions == pytest.approx(positions_before)


def test_crop_chains_agrees_with_apply_cropping(propyl):
    """
    The two entry points ran separate copies of the same capping loop and
    set different topology attributes -- `crop_chains` wrote `types` while
    `_apply_cropping` wrote `name`/`element`, so only one of the two
    returned a universe with a usable `atoms.elements`.
    """
    from_crop_chains, _keep, _replace = ChainCropper(
        cap_distance=CAP_DISTANCE).crop_chains(propyl, "alkyl", 1)

    processor = TrajectoryProcessor(cap_distance=CAP_DISTANCE)
    keep, delete, replace = processor.identify_chains_to_crop(propyl, "alkyl", 1)
    processor.keep_indices, processor.delete_indices, processor.replace_indices = \
        keep, delete, replace
    from_apply_cropping = processor._apply_cropping(propyl)

    assert (from_crop_chains.atoms.n_atoms
            == from_apply_cropping.atoms.n_atoms)
    assert list(from_crop_chains.atoms.names) == list(
        from_apply_cropping.atoms.names)
    assert list(from_crop_chains.atoms.elements) == list(
        from_apply_cropping.atoms.elements)
    assert from_crop_chains.atoms.positions == pytest.approx(
        from_apply_cropping.atoms.positions)


# --------------------------------------------------------------------------
# Trajectory paths
# --------------------------------------------------------------------------

def test_process_frame_agrees_with_apply_cropping(propyl):
    """
    `_process_frame` is the parallel trajectory path and `_apply_cropping`
    the serial one. They carried separate copies of the capping loop, which
    is how `_process_frame` kept a quadratic list-membership test long after
    it was fixed in `_apply_cropping`, and how neither of them got the
    one-cap-per-severed-bond fix. They must produce identical frames.
    """
    processor = TrajectoryProcessor(cap_distance=CAP_DISTANCE)
    keep, delete, replace = processor.identify_chains_to_crop(propyl, "alkyl", 1)
    processor.keep_indices, processor.delete_indices, processor.replace_indices = \
        keep, delete, replace

    serial = processor._apply_cropping(propyl)
    coords, elements, names = processor._process_frame(
        0,
        propyl.atoms.positions.copy(),
        np.asarray(propyl.atoms.types, dtype=object),
        np.asarray(propyl.atoms.names, dtype=object),
    )

    assert coords == pytest.approx(serial.atoms.positions)
    assert list(names) == list(serial.atoms.names)
    assert list(elements) == list(serial.atoms.elements)


def test_process_trajectory_writes_structure(propyl, tmp_path):
    output = tmp_path / "cropped.pdb"
    TrajectoryProcessor(cap_distance=CAP_DISTANCE).process_trajectory(
        structure_file=None, output_path=str(output),
        chain_type="alkyl", max_chain_length=1,
        structure_universe=propyl,
    )

    assert output.exists()
    written = mda.Universe(str(output))
    assert written.atoms.n_atoms == 5
