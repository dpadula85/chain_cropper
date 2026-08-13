#!/usr/bin/env python3
"""
Connectivity and hybridisation analysis on an MDAnalysis Universe.

This is the single home for the bond -> neighbour-matrix -> sp2/sp3
classification chain used across this pipeline. `SelIntCoords` and
`oligomer_builder` used to each carry their own copy of `get_sp2`; both now
re-export the implementation here, so a fix lands once.

The classification is deliberately connectivity-only -- no geometry, no
bond orders, no element-specific valence table -- because that is all the
upstream MD topologies reliably provide:

- **sp3** ("saturated"): a heavy atom with four or more connections. Such
  an atom cannot be part of a pi system, so it is what a side-chain walk
  is allowed to traverse.
- **sp2** ("unsaturated"): a heavy atom with fewer than four connections
  that is not itself a side-chain atom. In a cropped conjugated polymer
  this is the pi-conjugated core.

An atom with more than four connections (a sulfone sulfur, say) counts as
saturated, which is the same answer the old fixed-width-4 code gave -- but
it gets there without silently discarding the fifth bond.
"""

import numpy as np


def bond_indices(universe):
    """
    Bonds of `universe` as an (n_bonds, 2) integer array.

    Guesses bonds from geometry if the topology carries none. Guessing is
    both unreliable on a strained MD snapshot and the slow step on a large
    system, so prefer handing in a Universe built from a real topology.

    Parameters
    ----------
    universe : MDAnalysis.Universe

    Returns
    -------
    np.ndarray, shape (n_bonds, 2)
    """
    try:
        bonds = universe.bonds.to_indices()
    except (AttributeError, ValueError):
        # MDAnalysis raises NoDataError (a subclass of both) when the
        # topology has no bond information at all.
        universe.atoms.guess_bonds()
        bonds = universe.bonds.to_indices()

    return np.asarray(bonds, dtype=int)


def build_connectivity(universe, bonds=None, min_width=4):
    """
    Build a -1-padded neighbour matrix and the per-atom connection count.

    Parameters
    ----------
    universe : MDAnalysis.Universe
    bonds : np.ndarray, shape (n_bonds, 2), optional
        Bond list to use. Read (or guessed) off `universe` when omitted.
    min_width : int, default=4
        Minimum number of neighbour columns. The matrix is widened beyond
        this if some atom has more connections, so no bond is ever
        dropped. The old implementation hard-coded a width of 4 and
        silently discarded any further neighbour, which corrupted the
        chain walk and the capping around hypervalent atoms.

    Returns
    -------
    connectivity : np.ndarray, shape (n_atoms, width)
        Neighbour indices per atom, left-packed, `-1` in unused columns.
        Neighbours of a given atom appear in ascending order.
    degree : np.ndarray, shape (n_atoms,)
        Number of connections per atom. Read this rather than counting
        non-negative entries of a fixed-width row.
    """
    n_atoms = len(universe.atoms)

    if bonds is None:
        bonds = bond_indices(universe)
    else:
        bonds = np.asarray(bonds, dtype=int)

    if bonds.size == 0:
        return (np.full((n_atoms, min_width), -1, dtype=int),
                np.zeros(n_atoms, dtype=int))

    bonds = bonds.reshape(-1, 2)

    # Drop self-bonds and collapse duplicates: MDAnalysis topologies can
    # carry the same bond twice, which would otherwise consume two
    # neighbour slots and inflate the degree.
    bonds = bonds[bonds[:, 0] != bonds[:, 1]]
    bonds = np.unique(np.sort(bonds, axis=1), axis=0)

    # Both directions, grouped by owning atom.
    pairs = np.concatenate([bonds, bonds[:, ::-1]])
    pairs = pairs[np.lexsort((pairs[:, 1], pairs[:, 0]))]

    degree = np.bincount(pairs[:, 0], minlength=n_atoms)
    width = max(min_width, int(degree.max()))

    connectivity = np.full((n_atoms, width), -1, dtype=int)
    # Column of each neighbour within its owner's row: position in the
    # sorted list minus where that owner's block starts.
    starts = np.cumsum(degree) - degree
    column = np.arange(len(pairs)) - np.repeat(starts, degree)
    connectivity[pairs[:, 0], column] = pairs[:, 1]

    return connectivity, degree


def sp2_sp3(universe, alkyl=True, ether=False, connectivity=None, degree=None):
    """
    Split the heavy atoms into unsaturated (sp2) and saturated (sp3) sets.

    Parameters
    ----------
    universe : MDAnalysis.Universe
    alkyl : bool, default=True
        Kept for signature compatibility with the copies this replaces.
        Side chains are walked over saturated atoms whether or not this is
        set, so it has no effect on its own.
    ether : bool, default=False
        Also treat oxygen as a side-chain atom, so that an ether oxygen
        joining two saturated stretches does not read as part of the
        conjugated core. **In the copies this replaces, this argument had
        no effect**: they tested `if alkyl: ... elif ether: ...`, and
        `alkyl` defaults to True, so the `ether` branch was unreachable
        unless a caller also passed `alkyl=False`. Nobody did.
    connectivity, degree : np.ndarray, optional
        Output of `build_connectivity`, to avoid recomputing it.

    Returns
    -------
    sp2 : np.ndarray
        Indices of unsaturated heavy atoms that are not side-chain atoms.
    sp3 : np.ndarray
        Indices of saturated atoms (four or more connections).
    """
    if connectivity is None or degree is None:
        connectivity, degree = build_connectivity(universe)

    types = np.asarray(universe.atoms.types, dtype=str)

    heavy = np.flatnonzero(types != "H")
    sp3 = np.flatnonzero(degree >= 4)

    if ether:
        side_chain = np.union1d(sp3, np.flatnonzero(types == "O"))
    else:
        side_chain = sp3

    # Unsaturated heavy atoms, minus anything that counts as side chain.
    # The saturated atoms are already excluded by the test itself, so the
    # subtraction only bites in ether mode, where a two-coordinate oxygen
    # would otherwise be read as part of the pi system.
    unsaturated = np.flatnonzero(degree < 4)
    sp2 = np.intersect1d(np.setdiff1d(unsaturated, side_chain), heavy)

    return sp2, sp3


def get_sp2(u, alkyl=True, ether=False):
    """
    Backwards-compatible alias for `sp2_sp3`.

    Kept under this name because `oligomer_builder.enhanced_breaker`,
    `SelIntCoords.sel_intcoords` and `pyscf_tints.aom.overlap` all import
    a symbol called `get_sp2`.
    """
    return sp2_sp3(u, alkyl=alkyl, ether=ether)


def side_chain_atoms(universe, chain_type="alkyl", connectivity=None, degree=None):
    """
    Indices of the atoms a side-chain walk of the given type may traverse.

    A side chain is a *pendant* group: it hangs off the core at a single
    point. So a saturated atom is only side chain if at most one of its
    heavy neighbours belongs to the core. A saturated atom with two or more
    core neighbours is a **bridge** -- it holds two pieces of core together,
    and is therefore core itself.

    That distinction is what keeps fluorene, cyclopentadithiophene and
    indacenodithiophene intact. Their bridging carbons are four-coordinate,
    so a walk that treats every four-coordinate atom as side chain reaches
    them, and then either deletes them (splitting the conjugated unit into
    two uncapped halves) or leaves them classified as neither kept nor
    deleted, in which case they vanish from the output while the hydrogens
    that hung off them stay behind, bonded to nothing.

    An all-saturated ring hanging off the core by one bond -- a cyclohexyl
    substituent -- is still side chain, correctly: its attachment atom has
    exactly one core neighbour.

    Parameters
    ----------
    universe : MDAnalysis.Universe
    chain_type : {'alkyl', 'ether'}
        'alkyl' walks saturated atoms only. 'ether' additionally walks
        oxygen, so a glycol/alkoxy chain is followed through its oxygens
        instead of being chopped into disconnected saturated stretches.
        Note that an aryl-O-aryl oxygen is a bridge by the rule above, so
        turning ether mode on does not start cutting diaryl ethers in half.
    connectivity, degree : np.ndarray, optional
        Output of `build_connectivity`, to avoid recomputing it.

    Returns
    -------
    np.ndarray
        Sorted atom indices.

    Raises
    ------
    ValueError
        If `chain_type` is not recognised. The implementations this
        replaces silently fell back to 'alkyl' for any unknown value.
    """
    if chain_type not in ("alkyl", "ether"):
        raise ValueError(
            f"unknown chain_type {chain_type!r}; expected 'alkyl' or 'ether'"
        )

    if connectivity is None or degree is None:
        connectivity, degree = build_connectivity(universe)

    types = np.asarray(universe.atoms.types, dtype=str)
    n_atoms = len(types)
    is_heavy = types != "H"

    _sp2, sp3 = sp2_sp3(
        universe, ether=(chain_type == "ether"),
        connectivity=connectivity, degree=degree,
    )

    candidate = np.zeros(n_atoms, dtype=bool)
    candidate[sp3] = True
    if chain_type == "ether":
        candidate[types == "O"] = True

    # Core = every heavy atom that is not even a candidate side-chain atom.
    is_core = is_heavy & ~candidate

    valid = connectivity >= 0
    safe = np.where(valid, connectivity, 0)
    n_core_neighbours = (valid & is_core[safe]).sum(axis=1)

    return np.flatnonzero(candidate & (n_core_neighbours < 2))
