from scipy.spatial import cKDTree
import os
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdmolops
from ..rdutil import position
from ..rdutil import get_connected_component_indices
from ..rdutil import keep_atoms
from ..rdutil.select import match_smarts
from ..rdutil.geometry import check_occlusion
from ..rdutil import read_molecules
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

def write_orca_input(
    mol,
    center_atom_index,
    header,
    out_folder,
    mdci_header=None,
    additional_indices=None,
):
    if mdci_header is None:
        mdci_header = ""
    if additional_indices is None:
        additional_indices = []
    if header is None:
        header = "ORCA_HEADER_PLACEHOLDER"
    if type(mol) == str:
        mol = next(read_molecules(mol))
    cutout_mol, split_fragments, split_covalent_bonds = fragment_molecule(
        mol,
        center_atom_index,
        verbose=False,
        filter_occluded_fragments=False,
        additional_indices=additional_indices,
    )
    os.makedirs(out_folder, exist_ok=True)
    os.makedirs(os.path.join(out_folder, "super"), exist_ok=True)
    os.makedirs(os.path.join(out_folder, "sub1"), exist_ok=True)
    os.makedirs(os.path.join(out_folder, "sub2"), exist_ok=True)

    super_path = os.path.join(out_folder, "super", "orca.inp")
    sub1_path = os.path.join(out_folder, "sub1", "orca.inp")
    sub2_path = os.path.join(out_folder, "sub2", "orca.inp")

    write_fragment_input(
        cutout_mol,
        split_fragments,
        split_covalent_bonds,
        header,
        mdci_header,
        super_path,
        center_fragment_last=True,
    )

    if split_fragments:
        center_frag_index = len(split_fragments) - 1
        sub1_mol, sub1_frags, sub1_bonds = _subset_fragments(
            cutout_mol,
            split_fragments,
            split_covalent_bonds,
            [i for i in range(len(split_fragments)) if i != center_frag_index],
        )
        write_fragment_input(
            sub1_mol,
            sub1_frags,
            sub1_bonds,
            header,
            mdci_header,
            sub1_path,
            center_fragment_last=False,
        )

        sub2_mol, sub2_frags, sub2_bonds = _subset_fragments(
            cutout_mol,
            split_fragments,
            split_covalent_bonds,
            [center_frag_index],
        )
        write_fragment_input(
            sub2_mol,
            sub2_frags,
            sub2_bonds,
            header,
            mdci_header,
            sub2_path,
            center_fragment_last=False,
        )

def write_fragment_input(
    mol,
    fragment_indices,
    covalent_bonds,
    header,
    mdci_header,
    out_name,
    center_fragment_last=False,
):
    fragment_map = {}
    for frag_idx, frag in enumerate(fragment_indices, start=1):
        for atom_idx in frag:
            fragment_map[int(atom_idx)] = frag_idx
    ordered_bonds = _order_covalent_bonds(
        covalent_bonds, fragment_indices, center_fragment_last=center_fragment_last
    )
    lines = []
    lines.append(header)
    lines.extend(_format_mdci_block(ordered_bonds, mdci_header))
    lines.append("")
    lines.append("*xyz 0 1")
    positions = position(mol)
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        frag_idx = fragment_map.get(idx, 1)
        pos = positions[idx]
        lines.append(
            f"\t{atom.GetSymbol()}({frag_idx})\t{pos[0]:12.6f}\t{pos[1]:12.6f}\t{pos[2]:12.6f}"
        )
    lines.append("*")
    with open(out_name, "w") as f:
        f.write("\n".join(lines))

def _format_mdci_block(covalent_bonds, mdci_header):
    lines = ["%mdci"]
    if covalent_bonds:
        covalent_line = "  Covalent " + " ".join(
            f"{{{a} {b}}}" for a, b in covalent_bonds
        )
        lines.append(covalent_line)
    if mdci_header:
        for line in str(mdci_header).splitlines():
            if line.strip():
                lines.append(f"  {line}")
    lines.append("end")
    return lines

def _order_covalent_bonds(covalent_bonds, fragments, center_fragment_last):
    if not covalent_bonds:
        return []
    if not center_fragment_last or not fragments:
        return sorted(covalent_bonds)
    center_set = set(fragments[-1])
    first = []
    rest = []
    for a, b in covalent_bonds:
        if (a in center_set) ^ (b in center_set):
            first.append((a, b))
        else:
            rest.append((a, b))
    return sorted(first) + sorted(rest)

def _subset_fragments(mol, fragments, covalent_bonds, fragment_indices):
    if not fragment_indices:
        return Chem.Mol(), [], []
    selected_fragments = [fragments[i] for i in fragment_indices]
    keep_atoms_set = set()
    for frag in selected_fragments:
        keep_atoms_set.update(frag)
    keep_atoms_sorted = sorted(keep_atoms_set)
    index_map = {old_idx: new_idx for new_idx, old_idx in enumerate(keep_atoms_sorted)}
    new_mol = keep_atoms(mol, keep_atoms_sorted)
    new_fragments = []
    for frag in selected_fragments:
        new_fragments.append([index_map[idx] for idx in frag if idx in index_map])
    new_bonds = []
    for a, b in covalent_bonds:
        if a in index_map and b in index_map:
            new_bonds.append(tuple(sorted((index_map[a], index_map[b]))))
    return new_mol, new_fragments, sorted(new_bonds)

def fragment_molecule(
    mol,
    center_atom_index,
    verbose=False,
    filter_occluded_fragments=True,
    additional_indices=None,
):
    """
    Assumes the molecule is already proximity bonded.
    """
    all_components = get_connected_component_indices(mol)
    if verbose:
        logger.info(f"Connected components: {len(all_components)}")
    center_component_indices = None
    for component_indices in all_components:
        if center_atom_index in component_indices:
            center_component_indices = component_indices
            break
    if verbose:
        logger.info(f"Found center component {center_component_indices}")
    if center_component_indices is None:
        raise ValueError(f"Center atom index {center_atom_index} not found in molecule")

    surrounding_indices = get_surrounding_indices(mol, center_component_indices)
    if additional_indices:
        surrounding_indices = list(set(surrounding_indices).union(set(additional_indices)))
    if verbose:
        logger.info(f"Surrounding indices: {len(surrounding_indices)}")

    # extend surrounding indices to next C-C bond.
    all_used_components = []
    all_to_replace = []
    ring_info = None
    try:
        rdmolops.FastFindRings(mol)
        ring_info = mol.GetRingInfo()
    except Exception:
        ring_info = None

    for component_indices in all_components:
        current_selected = set(component_indices).intersection(surrounding_indices)
        current_selected, to_replace = extend_connected_indices(
            mol, current_selected, component_indices, ring_info=ring_info, max_ring_size=12
        )
        if len(current_selected) != 0:
            all_used_components.append(current_selected)
            if verbose:
                logger.info(f"Extended component {len(component_indices)} to {len(current_selected)}")
            all_to_replace.extend(to_replace)
    if verbose:
        logger.info(f"Extended components total: {len(all_used_components)}")

    split_fragments = split_fragments_by_peptide_bonds(all_used_components, mol, verbose=verbose)
    if verbose:
        logger.info(f"Split fragments: {len(split_fragments)}")

    if filter_occluded_fragments:
        split_fragments = _filter_fragments_by_occlusion(
            mol,
            split_fragments,
            center_atom_index,
            additional_indices=additional_indices,
            verbose=verbose,
        )
        if verbose:
            logger.info(f"Occlusion filter: {len(split_fragments)}")

    kept_atoms = set()
    for frag in split_fragments:
        kept_atoms.update(frag)
    to_replace = _compute_replacement_pairs(mol, split_fragments, kept_atoms)
    if verbose:
        logger.info(f"Replacement pairs: {len(to_replace)}")

    cutout_mol, index_map, added_hydrogens = create_molecule_cutout(
        mol, split_fragments, to_replace, verbose=verbose
    )
    if verbose:
        logger.info(f"Cutout mol atoms: {cutout_mol.GetNumAtoms()}")
    split_covalent_bonds = _get_fragment_connection_bonds(mol, split_fragments, index_map)
    split_fragments = [
        sorted(index_map[idx] for idx in frag if idx in index_map)
        for frag in split_fragments
        if frag
    ]
    if added_hydrogens:
        atom_to_fragment = {
            atom_idx: frag_idx
            for frag_idx, frag in enumerate(split_fragments)
            for atom_idx in frag
        }
        for base_old_idx, h_new_idx in added_hydrogens:
            base_new_idx = index_map.get(base_old_idx)
            if base_new_idx is None:
                continue
            frag_idx = atom_to_fragment.get(base_new_idx)
            if frag_idx is None:
                continue
            split_fragments[frag_idx].append(h_new_idx)

    center_new_idx = index_map.get(center_atom_index)
    cutout_mol, split_fragments, split_covalent_bonds = _reorder_fragments_consecutive(
        cutout_mol, split_fragments, split_covalent_bonds, center_new_idx
    )

    return cutout_mol, split_fragments, split_covalent_bonds

def _compute_replacement_pairs(mol, frags, kept_atoms):
    to_replace = set()
    for frag in frags:
        frag_set = set(frag)
        for idx in frag_set:
            atom = mol.GetAtomWithIdx(int(idx))
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in kept_atoms and neighbor.GetAtomicNum() != 1:
                    to_replace.add((idx, neighbor_idx))
    return list(to_replace)

def _get_fragment_connection_bonds(mol, frags, index_map):
    atom_to_fragment = {}
    for frag_idx, frag in enumerate(frags):
        for atom_idx in frag:
            atom_to_fragment[int(atom_idx)] = frag_idx
    connection_bonds = set()
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if a1 not in atom_to_fragment or a2 not in atom_to_fragment:
            continue
        if atom_to_fragment[a1] == atom_to_fragment[a2]:
            continue
        if a1 in index_map and a2 in index_map:
            connection_bonds.add(tuple(sorted((index_map[a1], index_map[a2]))))
    return sorted(connection_bonds)

def _reorder_fragments_consecutive(mol, fragments, covalent_bonds, center_atom_index):
    ordered_fragments = [sorted(frag) for frag in fragments]
    center_idx = None
    if center_atom_index is not None:
        for i, frag in enumerate(ordered_fragments):
            if center_atom_index in frag:
                center_idx = i
                break
    center_frag = None
    if center_idx is not None:
        center_frag = ordered_fragments.pop(center_idx)
    new_order = []
    for frag in ordered_fragments:
        new_order.extend(frag)
    total_atoms = mol.GetNumAtoms()
    if len(new_order) != total_atoms:
        center_set = set(center_frag or [])
        new_order_set = set(new_order)
        remaining = [
            idx for idx in range(total_atoms)
            if idx not in new_order_set and idx not in center_set
        ]
        if remaining:
            ordered_fragments.append(remaining)
            new_order.extend(remaining)
    if center_frag is not None:
        ordered_fragments.append(center_frag)
        new_order.extend(center_frag)
    index_map = {old_idx: new_idx for new_idx, old_idx in enumerate(new_order)}
    reordered_mol = Chem.RenumberAtoms(mol, new_order)
    new_fragments = []
    offset = 0
    for frag in ordered_fragments:
        size = len(frag)
        new_fragments.append(list(range(offset, offset + size)))
        offset += size
    new_covalent_bonds = sorted(
        tuple(sorted((index_map[a], index_map[b]))) for a, b in covalent_bonds
    )
    return reordered_mol, new_fragments, new_covalent_bonds

def _filter_fragments_by_occlusion(
    mol, fragments, center_atom_index, additional_indices=None, verbose=False
):
    central_frag = None
    for frag in fragments:
        if center_atom_index in frag:
            central_frag = frag
            break
    if central_frag is None:
        if verbose:
            logger.info("No central fragment found for occlusion filter")
        return fragments

    central_atoms = [mol.GetAtomWithIdx(int(i)) for i in central_frag]
    central_positions = position(central_atoms)
    tree = cKDTree(central_positions)
    all_frag_atoms = set()
    for frag in fragments:
        all_frag_atoms.update(frag)

    additional_set = set(additional_indices or [])
    kept = []
    removed = 0
    for frag in fragments:
        if frag is central_frag:
            kept.append(frag)
            continue
        if additional_set.intersection(frag):
            kept.append(frag)
            continue
        frag_atoms = [mol.GetAtomWithIdx(int(i)) for i in frag]
        frag_positions = position(frag_atoms)
        _, nearest_idxs = tree.query(frag_positions)
        to_atoms = [central_atoms[int(i)] for i in nearest_idxs]
        occluder_indices = all_frag_atoms - set(frag) - set(central_frag)
        if not occluder_indices:
            kept.append(frag)
            continue
        occluder_atoms = [mol.GetAtomWithIdx(int(i)) for i in occluder_indices]
        occluded = check_occlusion(frag_atoms, to_atoms, occluder_atoms, ignore_outside=True)
        if np.all(occluded):
            removed += 1
            if verbose:
                logger.info(f"Occlusion removed fragment of size {len(frag)}")
            continue
        kept.append(frag)
    if verbose:
        logger.info(f"Occlusion filter removed {removed} fragments")
    return kept

def create_molecule_cutout(mol, frags, to_replace, verbose=False):
    all_indices = set()
    for frag in frags:
        all_indices.update(frag)
    remaining_mol = Chem.RWMol(keep_atoms(mol, all_indices))
    if verbose:
        logger.info(
            f"keep_atoms: kept {len(all_indices)} / {mol.GetNumAtoms()} atoms"
        )
    index_map = {old_idx: new_idx for new_idx, old_idx in enumerate(sorted(list(all_indices)))}
    # replace to_replace atoms with hydrogens at distance 1.1 and bond order 1
    added_hydrogens = []
    for atom_pair in to_replace:
        atom1 = mol.GetAtomWithIdx(atom_pair[0])
        atom2 = mol.GetAtomWithIdx(atom_pair[1])
        remaining_mol.AddAtom(Chem.Atom(1))
        new_idx = remaining_mol.GetNumAtoms() - 1
        remaining_mol.AddBond(index_map[atom1.GetIdx()], new_idx, Chem.rdchem.BondType.SINGLE)
        pos1 = position(atom1)
        pos2 = position(atom2)
        direction = pos2 - pos1
        direction /= np.linalg.norm(direction)
        new_position = pos1 + direction * 1.1
        remaining_mol.GetConformer().SetAtomPosition(new_idx, new_position)
        added_hydrogens.append((atom1.GetIdx(), new_idx))
    if verbose:
        logger.info(f"Hydrogen caps: {len(to_replace)}")
    return remaining_mol, index_map, added_hydrogens

def get_surrounding_indices(mol, center_indices, radius = 5):
    positions = position(mol)
    tree = cKDTree(positions)
    surrounding_indices = tree.query_ball_point(positions[center_indices], radius)
    surrounding_indices = [index for sublist in surrounding_indices for index in sublist]
    return surrounding_indices

def extend_connected_indices(mol, indices, component_indices, ring_info=None, max_ring_size=12):
    """ Extends the indices within the connected component of the given indices to the next C-C bond."""
    to_check = list(indices)
    seen = set(indices)
    keep = set(indices)
    component_set = set(component_indices)
    small_ring_atoms = set()
    if ring_info is not None:
        try:
            for idx in component_set:
                for size in range(3, max_ring_size + 1):
                    if ring_info.IsAtomInRingOfSize(int(idx), int(size)):
                        small_ring_atoms.add(idx)
                        break
        except Exception:
            small_ring_atoms = set()
    # heavy atom check
    while to_check:
        index = to_check.pop()
        index_atom = mol.GetAtomWithIdx(index)
        for neighbor in mol.GetAtomWithIdx(index).GetNeighbors():
            neighbor_index = neighbor.GetIdx()
            if neighbor_index not in component_set:
                continue
            neighbor_symbol = neighbor.GetSymbol()
            index_symbol = index_atom.GetSymbol()
            single_bond = mol.GetBondBetweenAtoms(index, neighbor_index).GetBondType() == Chem.rdchem.BondType.SINGLE
            both_carbon = (index_symbol == "C") and (neighbor_symbol == "C")
            # do not cut rings
            if not (neighbor_index in seen) and (not(single_bond and both_carbon and not(neighbor_index in small_ring_atoms) and not(index in small_ring_atoms))):
                to_check.append(neighbor_index)
            seen.add(neighbor_index)
            if both_carbon:
                # check if neighbor is sp3, we proxy by having 4 neighbors since we know its carbon
                if len(neighbor.GetNeighbors()) == 4:
                    keep.add(neighbor_index)
            else:
                keep.add(neighbor_index)
    # add any hydrogen neighbors
    new_seen = set()
    for index in keep:
        for neighbor in mol.GetAtomWithIdx(index).GetNeighbors():
            neighbor_index = neighbor.GetIdx()
            if neighbor.GetAtomicNum() == 1:
                new_seen.add(neighbor_index)
    
    # determine atoms to replace with hydrogens
    to_replace = set()
    for index in keep:
        for neighbor in mol.GetAtomWithIdx(index).GetNeighbors():
            neighbor_index = neighbor.GetIdx()
            if neighbor_index not in keep:
                if neighbor.GetAtomicNum() != 1:
                    to_replace.add((index, neighbor_index))
    keep.update(new_seen)
    return list(keep), list(to_replace)

def _get_mapped_atom_indices(pattern, map_numbers):
    mapped = {}
    for atom in pattern.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num in map_numbers:
            mapped[map_num] = atom.GetIdx()
    missing = [num for num in map_numbers if num not in mapped]
    if missing:
        raise ValueError(f"SMARTS pattern missing atom map numbers: {missing}")
    return mapped

def _find_bond_pairs_from_smarts(mol, smarts, map_num_a=1, map_num_b=2):
    matches = match_smarts(mol, smarts)
    if not matches:
        return []
    bond_pairs = set()
    for match in matches:
        atom_a = match[map_num_a - 1]
        atom_b = match[map_num_b - 1]
        bond_pairs.add(tuple(sorted((atom_a, atom_b))))
    return sorted(bond_pairs)

def _split_fragment_on_bond(fragment_set, bond_pair, fragment_mask, neighbor_lists):
    atom_a, atom_b = bond_pair
    if atom_a not in fragment_set or atom_b not in fragment_set:
        return None
    visited = set()
    stack = [atom_a]
    while stack:
        current = stack.pop()
        if current in visited:
            continue
        visited.add(current)
        for neighbor in neighbor_lists.get(current, []):
            if neighbor not in fragment_mask:
                continue
            if (current == atom_a and neighbor == atom_b) or (current == atom_b and neighbor == atom_a):
                continue
            if neighbor not in visited:
                stack.append(neighbor)
    if atom_b in visited:
        return None
    other = fragment_set - visited
    if not other:
        return None
    return visited, other

def _iterative_split_by_bonds(fragment_set, bond_pairs_by_atom, min_size, neighbor_lists):
    pending = [fragment_set]
    final_fragments = []
    while pending:
        current = pending.pop()
        fragment_mask = set(current)
        split_done = False
        candidate_pairs = set()
        for atom in current:
            for pair in bond_pairs_by_atom.get(atom, []):
                candidate_pairs.add(pair)
        for bond_pair in candidate_pairs:
            if bond_pair[0] not in current or bond_pair[1] not in current:
                continue
            split = _split_fragment_on_bond(current, bond_pair, fragment_mask, neighbor_lists)
            if split is None:
                continue
            left, right = split
            if len(left) > min_size and len(right) > min_size:
                pending.append(left)
                pending.append(right)
                split_done = True
                break
        if not split_done:
            final_fragments.append(current)
    return final_fragments

def _build_bond_pair_map(bond_pairs):
    bond_pairs_by_atom = {}
    for pair in bond_pairs:
        atom_a, atom_b = sorted(pair)
        bond_pairs_by_atom.setdefault(atom_a, []).append((atom_a, atom_b))
        bond_pairs_by_atom.setdefault(atom_b, []).append((atom_a, atom_b))
    return bond_pairs_by_atom

def _split_into_connected_components(mol, fragment_set):
    components = []
    remaining = set(fragment_set)
    while remaining:
        start = next(iter(remaining))
        stack = [start]
        visited = set()
        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)
            for neighbor in mol.GetAtomWithIdx(int(current)).GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx in remaining and neighbor_idx not in visited:
                    stack.append(neighbor_idx)
        components.append(visited)
        remaining -= visited
    return components

def split_fragments_by_peptide_bonds(all_used_components, mol, verbose=False):
    """
    Further splits fragments by sidechain and backbone C-C bonds.

    Sidechain C-C bond: C-N-C(=O)-C (split at C(=O)-C).
    Backbone C-C bond:  C(=O)(N)-C(N) (split at C(=O)-C).
    """
    sidechain_smarts = "[C:1]([N]-[C](=O)-[C])-[C:2]"
    backbone_smarts = "[C:1](=O)([N])-[C:2]([N])"
    sidechain_bonds = _find_bond_pairs_from_smarts(mol, sidechain_smarts, 1, 2)
    backbone_bonds = _find_bond_pairs_from_smarts(mol, backbone_smarts, 1, 2)
    if verbose:
        logger.info(
            f"SMARTS matches: sidechain={len(sidechain_bonds)}, backbone={len(backbone_bonds)}"
        )
    sidechain_bonds_by_atom = _build_bond_pair_map(sidechain_bonds)
    backbone_bonds_by_atom = _build_bond_pair_map(backbone_bonds)
    all_used_atoms = set()
    for frag in all_used_components:
        all_used_atoms.update(frag)

    sidechain_bonds = [pair for pair in sidechain_bonds if (pair[0] in all_used_atoms and pair[1] in all_used_atoms)]
    backbone_bonds = [pair for pair in backbone_bonds if (pair[0] in all_used_atoms and pair[1] in all_used_atoms)]
    if verbose:
        logger.info(
            f"Filtered bonds: sidechain={len(sidechain_bonds)}, backbone={len(backbone_bonds)}"
        )
    neighbor_lists = {}
    for idx in all_used_atoms:
        neighbor_lists[idx] = [n.GetIdx() for n in mol.GetAtomWithIdx(int(idx)).GetNeighbors()]

    final_fragments = []
    for frag in all_used_components:
        fragment_set = set(frag)
        sidechain_split = _iterative_split_by_bonds(
            fragment_set, sidechain_bonds_by_atom, min_size=3, neighbor_lists=neighbor_lists
        )
        for subfrag in sidechain_split:
            backbone_split = _iterative_split_by_bonds(
                subfrag, backbone_bonds_by_atom, min_size=4, neighbor_lists=neighbor_lists
            )
            for backbone_frag in backbone_split:
                final_fragments.extend(
                    _split_into_connected_components(mol, backbone_frag)
                )
    return [sorted(list(fragment)) for fragment in final_fragments]
