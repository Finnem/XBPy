"""
LED (Local Energy Decomposition) extraction and analysis utilities.

This module provides functions for extracting and processing LED matrices from ORCA
calculations, implementing the CovaLED and fp-CovaLED methodologies for proper
interaction energy decomposition.

LEDAW Integration
-----------------
This module bundles LEDAW (LED Analysis Workflow) as a vendored dependency in
`xbpy._vendor.ledaw_package`. LEDAW is developed and maintained by its original
authors - all credit belongs to the LEDAW development team.

The bundled copy is included for convenience to avoid installation/path issues.
See `xbpy/_vendor/VENDOR_README.md` for more information.

**No additional setup required** - LEDAW is imported automatically from the
vendored location

CovaLED Methodology
-------------------
The CovaLED (Covalent LED) approach properly decomposes interaction energies by
accounting for electronic reorganization within subsystems. See:
https://www.faccts.de/docs/opi/nightly/docs/contents/notebooks/covaled_ethane.html

Raw LED matrices from LEDAW contain:
- Inter-molecular interactions (what we want)
- Intra-molecular electronic reorganization/deformation

Simply summing inter-molecular terms gives incorrect interaction energy!

CovaLED Step 7: Subtract Subsystem Contributions
-------------------------------------------------
For each fragment pair (i,j) in the supersystem:
  - If both in ligand: LED_int[i,j] = LED_super[i,j] - LED_ligand[i',j']
  - If both in receptor: LED_int[i,j] = LED_super[i,j] - LED_receptor[i',j']
  - If inter-molecular: LED_int[i,j] = LED_super[i,j] (pure interaction)

After Step 7: Sum(ALL LED_int terms) = E_int
  - Inter-molecular terms = pure ligand-receptor interactions
  - Intra-molecular terms = electronic reorganization within each molecule

fp-CovaLED Step 8: Redistribute Intra-molecular to Inter-molecular
-------------------------------------------------------------------
To express E_int as ONLY fragment-pairwise interactions:
1. Calculate total intra-molecular contributions (ligand + receptor)
2. Redistribute proportionally to inter-molecular pairs: weight = |value|/Σ|values|
3. New value = original + (weight × intra_total)

After Step 8: Sum(inter-molecular pairs ONLY) = E_int
  - All contributions are ligand-receptor fragment pairs
  - Properly accounts for electronic reorganization in the fragment contributions
"""

import pandas as pd
import numpy as np
import os
import sys
from typing import Dict, List, Tuple, Optional, Union
from pathlib import Path


__all__ = [
    'compute_led_interaction_matrix',
    'extract_interaction_energy',
    'parse_orca_fragments',
    'determine_fragment_groups',
    'clear_subsystem_cache',
    'compute_fp_led_interactions',
    'compute_fp_led_interactions_with_mapping'
]


# ==============================================================================
# LEDAW IMPORT (Vendored)
# ==============================================================================
# Import LEDAW from vendored location
# LEDAW is bundled in xbpy._vendor.ledaw_package for ease of use

try:
    from xbpy._vendor.ledaw_package import engine_LED_N_body as _LEDAW_ENGINE
    _LEDAW_AVAILABLE = True
except ImportError as e:
    _LEDAW_ENGINE = None
    _LEDAW_AVAILABLE = False
    _LEDAW_IMPORT_ERROR = str(e)


# ==============================================================================
# ORCA OUTPUT PARSING
# ==============================================================================

def parse_orca_fragments(orca_out_path: Union[str, Path]) -> Tuple[List[int], List[List[float]]]:
    """
    Parse ORCA output file to extract fragment assignments and atomic positions.
    
    Extracts the fragment assignments from the INPUT FILE section of the ORCA output.
    
    Example from ORCA output:
        | 11>   C(1)       -0.673900        0.537400        0.576600
        | 12>   S(2)       -0.162900        2.270400        0.534600
        | 13>   H(1)       -0.247100        0.010400       -0.276800
        
    This indicates:
        - Atom 0 (C) belongs to fragment 1 (index 0 in 0-based) at position (-0.6739, 0.5374, 0.5766)
        - Atom 1 (S) belongs to fragment 2 (index 1 in 0-based) at position (-0.1629, 2.2704, 0.5346)
        - Atom 2 (H) belongs to fragment 1 (index 0 in 0-based) at position (-0.2471, 0.0104, -0.2768)
    
    Args:
        orca_out_path: Path to ORCA .out file
        
    Returns:
        Tuple of (atom_to_fragment, positions) where:
        - atom_to_fragment: List where index is atom number and value is fragment index (0-based)
        - positions: List of [x, y, z] coordinates for each atom in Angstroms
        
    Raises:
        FileNotFoundError: If ORCA output file doesn't exist
        ValueError: If fragment information cannot be parsed
        
    Example:
        >>> atom_frags, positions = parse_orca_fragments("complex.out")
        >>> # atom_frags = [0, 0, 0, 1, 1, 1, 2, 2]  # Atoms 0-2 in frag 0, 3-5 in frag 1, 6-7 in frag 2
        >>> # positions = [[x1,y1,z1], [x2,y2,z2], ...]
    """
    orca_out_path = Path(orca_out_path)
    
    if not orca_out_path.exists():
        raise FileNotFoundError(f"ORCA output file not found: {orca_out_path}")
    
    atom_to_fragment = []
    positions = []
    in_input_section = False
    in_xyz_section = False
    
    with open(orca_out_path, 'r') as f:
        for line in f:
            # Look for the INPUT FILE section
            if "INPUT FILE" in line:
                in_input_section = True
                continue
            
            # Look for start of xyz coordinates
            if in_input_section and "*xyz" in line:
                in_xyz_section = True
                continue
            
            # Parse xyz section
            # Format: | 11>   C(1)       -0.673900        0.537400        0.576600
            if in_xyz_section and "|" in line:
                # Extract the actual line content after the pipe
                parts = line.split("|", 1)
                if len(parts) < 2:
                    continue
                
                content = parts[1].strip()
                if not content:
                    continue
                
                # Remove the line number prefix (e.g., "11> ")
                if ">" in content:
                    content = content.split(">", 1)[1].strip()
                
                # Check for end of xyz section (line with just "*")
                if content == "*":
                    break
                
                if not content:
                    continue
                
                # Parse atom specification like "C(1)" or "H(5)" with coordinates
                tokens = content.split()
                if len(tokens) >= 4:  # Need atom spec + x, y, z
                    atom_spec = tokens[0]
                    # Extract fragment number from parentheses
                    if "(" in atom_spec and ")" in atom_spec:
                        try:
                            start = atom_spec.index("(")
                            end = atom_spec.index(")")
                            frag_num = int(atom_spec[start+1:end])
                            # Convert to 0-based indexing
                            atom_to_fragment.append(frag_num - 1)
                            
                            # Extract coordinates
                            x = float(tokens[1])
                            y = float(tokens[2])
                            z = float(tokens[3])
                            positions.append([x, y, z])
                        except (ValueError, IndexError):
                            # Not a valid atom line, skip
                            continue
    
    if not atom_to_fragment:
        raise ValueError(f"Could not parse fragment assignments from ORCA output: {orca_out_path}")
    
    return atom_to_fragment, positions


def determine_fragment_groups(
    atom_to_fragment: List[int],
    ligand_atom_indices: Optional[List[int]] = None,
    assign_mixed: str = 'error'
) -> Tuple[List[int], List[int]]:
    """
    Determine which fragments belong to ligand vs receptor based on atom assignments.
    
    Properly handles fragments that may contain atoms from both ligand and receptor.
    
    Args:
        atom_to_fragment: List mapping atom indices to fragment indices (0-based)
        ligand_atom_indices: List of atom indices that belong to ligand.
                            If None, assumes first N/2 atoms are ligand.
        assign_mixed: How to handle fragments with atoms from both ligand and receptor:
                     - 'error': Raise ValueError (default, safest)
                     - 'majority': Assign to whichever has more atoms
                     - 'ligand': Assign mixed fragments to ligand
                     - 'receptor': Assign mixed fragments to receptor
    
    Returns:
        Tuple of (ligand_fragments, receptor_fragments) where each is a sorted list
        of fragment indices in 1-based indexing (for LED matrix lookup)
        
    Raises:
        ValueError: If a fragment contains atoms from both ligand and receptor
                   and assign_mixed='error'
        
    Example:
        >>> atom_frags = [0, 0, 0, 1, 1, 1, 2, 2]  # 8 atoms in 3 fragments
        >>> ligand_atoms = [0, 1, 2, 3, 4, 5]      # First 6 atoms are ligand
        >>> lig_frags, rec_frags = determine_fragment_groups(atom_frags, ligand_atoms)
        >>> # lig_frags = [1, 2]  # Fragments 0,1 in 1-based = [1, 2]
        >>> # rec_frags = [3]     # Fragment 2 in 1-based = [3]
    """
    if ligand_atom_indices is None:
        # Default: assume first half of atoms are ligand
        ligand_atom_indices = list(range(len(atom_to_fragment) // 2))
    
    ligand_atom_set = set(ligand_atom_indices)
    receptor_atom_set = set(range(len(atom_to_fragment))) - ligand_atom_set
    
    # Group atoms by fragment
    fragment_atoms = {}
    for atom_idx, frag_idx in enumerate(atom_to_fragment):
        if frag_idx not in fragment_atoms:
            fragment_atoms[frag_idx] = []
        fragment_atoms[frag_idx].append(atom_idx)
    
    ligand_fragments = []
    receptor_fragments = []
    
    # Determine assignment for each fragment
    for frag_idx, atoms_in_frag in sorted(fragment_atoms.items()):
        atoms_in_frag_set = set(atoms_in_frag)
        ligand_atoms_in_frag = atoms_in_frag_set & ligand_atom_set
        receptor_atoms_in_frag = atoms_in_frag_set & receptor_atom_set
        
        has_ligand = len(ligand_atoms_in_frag) > 0
        has_receptor = len(receptor_atoms_in_frag) > 0
        
        # Convert to 1-based indexing for LED matrix
        frag_1based = frag_idx + 1
        
        if has_ligand and has_receptor:
            # Fragment contains atoms from both molecules
            if assign_mixed == 'error':
                raise ValueError(
                    f"Fragment {frag_1based} (0-based: {frag_idx}) contains atoms from both ligand and receptor!\n"
                    f"  Ligand atoms in fragment: {sorted(ligand_atoms_in_frag)}\n"
                    f"  Receptor atoms in fragment: {sorted(receptor_atoms_in_frag)}\n"
                    f"  This indicates the fragment definition spans both molecules.\n"
                    f"  You may need to use finer fragment definitions in ORCA input,\n"
                    f"  or set assign_mixed='majority'/'ligand'/'receptor' to force assignment."
                )
            elif assign_mixed == 'majority':
                if len(ligand_atoms_in_frag) >= len(receptor_atoms_in_frag):
                    ligand_fragments.append(frag_1based)
                else:
                    receptor_fragments.append(frag_1based)
            elif assign_mixed == 'ligand':
                ligand_fragments.append(frag_1based)
            elif assign_mixed == 'receptor':
                receptor_fragments.append(frag_1based)
            else:
                raise ValueError(f"Invalid assign_mixed value: {assign_mixed}")
        elif has_ligand:
            ligand_fragments.append(frag_1based)
        elif has_receptor:
            receptor_fragments.append(frag_1based)
    
    return sorted(ligand_fragments), sorted(receptor_fragments)


# ==============================================================================
# LEDAW INTEGRATION
# ==============================================================================

def _run_ledaw(
    orca_out_path: Union[str, Path],
    method: str = "DLPNO-CCSD(T)",
    conversion_factor: float = 2625.5,
    output_dir: Optional[Union[str, Path]] = None,
    verbose: bool = True
) -> Path:
    """
    Run LEDAW on an ORCA output file to generate LED matrices.
    
    Args:
        orca_out_path: Path to ORCA .out file
        method: QM method used ("DLPNO-CCSD(T)", "DLPNO-CCSD", "HFLD")
        conversion_factor: Energy conversion factor (default: 2625.5 for Hartree to kJ/mol)
        output_dir: Directory for LEDAW output (default: same as ORCA output)
        verbose: If True, print detailed progress information
    
    Returns:
        Path to the generated LED Excel file
        
    Raises:
        ImportError: If LEDAW package is not available
        FileNotFoundError: If ORCA output file doesn't exist
        RuntimeError: If LEDAW execution fails
    """
    orca_out_path = Path(orca_out_path)
    
    if verbose:
        print(f"\n{'='*80}")
        print(f"LEDAW EXECUTION - DEBUG MODE")
        print(f"{'='*80}")
        print(f"[CHECKPOINT 1] Starting LEDAW execution")
        print(f"  Input file: {orca_out_path}")
        print(f"  Method: {method}")
        print(f"  Conversion factor: {conversion_factor}")
    
    if not orca_out_path.exists():
        raise FileNotFoundError(f"ORCA output file not found: {orca_out_path}")
    
    if verbose:
        file_size = orca_out_path.stat().st_size
        print(f"[CHECKPOINT 2] File exists, size: {file_size / 1024:.2f} KB")
    
    if output_dir is None:
        output_dir = orca_out_path.parent
    else:
        output_dir = Path(output_dir)
    
    if verbose:
        print(f"[CHECKPOINT 3] Output directory: {output_dir}")
    
    # Check LEDAW availability
    if verbose:
        print(f"[CHECKPOINT 4] Checking LEDAW availability...")
    
    if not _LEDAW_AVAILABLE:
        error_msg = (
            "\n" + "="*80 + "\n"
            "ERROR: LEDAW not available\n"
            "="*80 + "\n"
            f"Vendored LEDAW import failed: {_LEDAW_IMPORT_ERROR}\n"
            "\n"
            "LEDAW is bundled in xbpy._vendor.ledaw_package but failed to import.\n"
            "This usually indicates:\n"
            "  - Incompatible compiled extensions for your Python/system\n"
            "  - Missing system libraries (BLAS, LAPACK, MKL, etc.)\n"
            "  - Corrupted vendored package files\n"
            "\n"
            "Try reinstalling xbpy or contact maintainers.\n"
            "="*80
        )
        raise ImportError(error_msg)
    
    if verbose:
        print(f"[CHECKPOINT 5] ✓ LEDAW is available (vendored)")
    
    # Prepare arguments
    if verbose:
        print(f"[CHECKPOINT 6] Preparing LEDAW arguments...")
        print(f"  main_filenames: {[str(orca_out_path)]}")
        print(f"  alternative_filenames: {['']}")
        print(f"  conversion_factor: {conversion_factor}")
        print(f"  method: {method}")
        print(f"  use_ref_as_rhf_in_hfld: {False}")
        print(f"  LEDAW_output_path: {str(output_dir)}")
    
    # Run LEDAW
    if verbose:
        print(f"[CHECKPOINT 7] Calling engine_LED_N_body (vendored LEDAW)...")
        print(f"{'='*80}")
        print(f"LEDAW OUTPUT STARTS BELOW:")
        print(f"{'='*80}")
    
    try:
        _LEDAW_ENGINE(
            main_filenames=[str(orca_out_path)],
            alternative_filenames=[""],
            conversion_factor=conversion_factor,
            method=method,
            use_ref_as_rhf_in_hfld=False,
            LEDAW_output_path=str(output_dir)
        )
        
        if verbose:
            print(f"{'='*80}")
            print(f"LEDAW OUTPUT ENDS ABOVE")
            print(f"{'='*80}")
            print(f"[CHECKPOINT 8] engine_LED_N_body returned successfully")
            
    except Exception as e:
        if verbose:
            print(f"{'='*80}")
            print(f"[CHECKPOINT 8 - ERROR] Exception caught!")
            print(f"  Exception type: {type(e).__name__}")
            print(f"  Exception message: {e}")
            import traceback
            traceback.print_exc()
            print(f"{'='*80}")
        raise RuntimeError(f"LEDAW execution failed for {orca_out_path}: {e}")
    
    # Check that Excel file was created
    if verbose:
        print(f"[CHECKPOINT 9] Checking for output Excel file...")
    
    excel_path = output_dir / "All_Standard_LED_matrices.xlsx"
    if not excel_path.exists():
        if verbose:
            print(f"[CHECKPOINT 9 - ERROR] Excel file not found!")
            print(f"  Expected: {excel_path}")
            print(f"  Directory contents:")
            for item in output_dir.iterdir():
                print(f"    - {item.name}")
        raise RuntimeError(f"LEDAW completed but Excel file not created: {excel_path}")
    
    if verbose:
        print(f"[CHECKPOINT 10] Excel file created successfully: {excel_path.name}")
        print(f"{'='*80}")
        print(f"LEDAW EXECUTION COMPLETED SUCCESSFULLY")
        print(f"{'='*80}\n")
    
    return excel_path


def _ensure_led_file(
    orca_out_path: Union[str, Path],
    method: str = "DLPNO-CCSD(T)",
    conversion_factor: float = 2625.5,
    force_regenerate: bool = False,
    verbose: bool = False
) -> Path:
    """
    Ensure LED Excel file exists for an ORCA output, running LEDAW if needed.
    
    Args:
        orca_out_path: Path to ORCA .out file
        method: QM method used
        conversion_factor: Energy conversion factor
        force_regenerate: If True, regenerate even if Excel file exists
        verbose: If True, enable detailed progress output
    
    Returns:
        Path to LED Excel file
    """
    orca_out_path = Path(orca_out_path)
    output_dir = orca_out_path.parent
    excel_path = output_dir / "All_Standard_LED_matrices.xlsx"
    
    if not excel_path.exists() or force_regenerate:
        if not verbose:
            print(f"  Generating LED matrices: {orca_out_path.name}")
        excel_path = _run_ledaw(orca_out_path, method, conversion_factor, output_dir, verbose=verbose)
        if not verbose:
            print(f"  ✓ LED matrices generated: {excel_path}")
    
    return excel_path


# ==============================================================================
# GLOBAL CACHE FOR SUBSYSTEM LED MATRICES
# ==============================================================================
# Since subsystem LED files are the same across many calculations, we cache them
# to avoid re-loading 100s of MB of Excel data repeatedly
_SUBSYSTEM_LED_CACHE: Dict[Tuple[str, str], pd.DataFrame] = {}


def _get_cached_subsystem_led(excel_path: Union[str, Path], component: str = 'TOTAL') -> pd.DataFrame:
    """
    Load and cache subsystem LED matrices to avoid repeated file reads.
    
    Args:
        excel_path: Path to LED Excel file
        component: LED component sheet name (default: 'TOTAL')
    
    Returns:
        DataFrame with LED matrix (cached after first load)
    """
    cache_key = (str(excel_path), component)
    
    if cache_key not in _SUBSYSTEM_LED_CACHE:
        _SUBSYSTEM_LED_CACHE[cache_key] = pd.read_excel(excel_path, sheet_name=component, index_col=0)
    
    return _SUBSYSTEM_LED_CACHE[cache_key]


def clear_subsystem_cache():
    """Clear the subsystem LED matrix cache to free memory."""
    global _SUBSYSTEM_LED_CACHE
    _SUBSYSTEM_LED_CACHE.clear()


def _compute_covaled_interaction_matrix(
    df_super: pd.DataFrame,
    df_ligand: pd.DataFrame,
    df_receptor: pd.DataFrame,
    ligand_fpled_frags: List[int],
    receptor_fpled_frags: List[int]
) -> pd.DataFrame:
    """
    Compute LED interaction matrix following CovaLED Step 7.
    
    Subtracts subsystem LED contributions from supersystem to get pure interaction energy:
    - Intra-ligand pairs: LED_int[i,j] = LED_super[i,j] - LED_ligand[i',j']
    - Intra-receptor pairs: LED_int[i,j] = LED_super[i,j] - LED_receptor[i',j']
    - Inter-molecular pairs: LED_int[i,j] = LED_super[i,j] (no subtraction)
    
    Args:
        df_super: Supersystem LED matrix (from complex calculation)
        df_ligand: Ligand subsystem LED matrix (from isolated ligand)
        df_receptor: Receptor subsystem LED matrix (from isolated receptor)
        ligand_fpled_frags: List of fragment indices (1-based) belonging to ligand
        receptor_fpled_frags: List of fragment indices (1-based) belonging to receptor
    
    Returns:
        DataFrame with LED interaction matrix (same indexing as supersystem)
        
    Example:
        >>> # Load LED matrices
        >>> df_super = pd.read_excel("complex_led.xlsx", sheet_name="TOTAL", index_col=0)
        >>> df_lig = pd.read_excel("ligand_led.xlsx", sheet_name="TOTAL", index_col=0)
        >>> df_rec = pd.read_excel("receptor_led.xlsx", sheet_name="TOTAL", index_col=0)
        >>> 
        >>> # Define fragment assignments (1-based indexing)
        >>> ligand_frags = [1, 2, 3]  # Fragments 1-3 belong to ligand
        >>> receptor_frags = [4, 5]   # Fragments 4-5 belong to receptor
        >>> 
        >>> # Compute interaction matrix
        >>> led_int = compute_covaled_interaction_matrix(
        ...     df_super, df_lig, df_rec, ligand_frags, receptor_frags
        ... )
        >>> 
        >>> # Sum all terms to get total interaction energy
        >>> e_int = led_int.sum().sum()
    """
    # Create interaction LED matrix (copy supersystem structure)
    led_int = df_super.copy()
    
    all_frags = sorted(set(ligand_fpled_frags + receptor_fpled_frags))
    
    # Process each fragment pair
    for i, frag_i in enumerate(all_frags):
        for j, frag_j in enumerate(all_frags):
            if frag_i > frag_j:  # Skip lower triangle
                continue
            
            try:
                super_val = df_super.loc[frag_i, frag_j]
                if pd.isna(super_val):
                    continue
                
                # Determine which molecule each fragment belongs to
                in_ligand_i = frag_i in ligand_fpled_frags
                in_ligand_j = frag_j in ligand_fpled_frags
                
                # Case 1: Both fragments in ligand - subtract ligand subsystem contribution
                if in_ligand_i and in_ligand_j:
                    # Map supersystem indices to subsystem indices (1-based)
                    subsys_i = ligand_fpled_frags.index(frag_i) + 1
                    subsys_j = ligand_fpled_frags.index(frag_j) + 1
                    try:
                        subsys_val = df_ligand.loc[subsys_i, subsys_j]
                        if pd.notna(subsys_val):
                            led_int.loc[frag_i, frag_j] = super_val - subsys_val
                    except KeyError:
                        # Subsystem value not found, keep supersystem value
                        pass
                
                # Case 2: Both fragments in receptor - subtract receptor subsystem contribution
                elif not in_ligand_i and not in_ligand_j:
                    subsys_i = receptor_fpled_frags.index(frag_i) + 1
                    subsys_j = receptor_fpled_frags.index(frag_j) + 1
                    try:
                        subsys_val = df_receptor.loc[subsys_i, subsys_j]
                        if pd.notna(subsys_val):
                            led_int.loc[frag_i, frag_j] = super_val - subsys_val
                    except KeyError:
                        pass
                
                # Case 3: Inter-molecular (one in ligand, one in receptor)
                # Keep supersystem value as-is (pure interaction, no subtraction needed)
                else:
                    pass  # led_int already has super_val from copy
                    
            except KeyError:
                continue
    
    return led_int


def _compute_fpcovaled_interactions(
    led_int: pd.DataFrame,
    ligand_fpled_frags: List[int],
    receptor_fpled_frags: List[int],
    fragment_pair_values: Dict[str, float]
) -> Dict[str, float]:
    """
    Apply fp-CovaLED Step 8: redistribute intra-molecular contributions to inter-molecular pairs.
    
    After CovaLED Step 7, the LED_int matrix contains:
    - Inter-molecular terms: pure ligand-receptor interactions
    - Intra-molecular terms: electronic reorganization within each molecule
    
    Sum(all LED_int) = E_int, but only inter-molecular pairs represent physical interactions.
    
    This function redistributes the intra-molecular contributions proportionally to
    inter-molecular pairs based on their absolute magnitudes, so that:
    Sum(inter-molecular pairs ONLY) = E_int
    
    Args:
        led_int: LED interaction matrix from CovaLED Step 7
        ligand_fpled_frags: List of fragment indices (1-based) belonging to ligand
        receptor_fpled_frags: List of fragment indices (1-based) belonging to receptor
        fragment_pair_values: Dict of inter-molecular fragment pair interactions
                             (before fp-CovaLED redistribution)
                             Format: {'frag_i_frag_j': energy, ...}
    
    Returns:
        Dict with fp-CovaLED redistributed interactions (same format as input)
        
    Example:
        >>> # After CovaLED Step 7
        >>> led_int = compute_covaled_interaction_matrix(...)
        >>> 
        >>> # Extract inter-molecular pairs (your application-specific logic)
        >>> inter_pairs = extract_intermolecular_pairs(led_int, lig_frags, rec_frags)
        >>> # inter_pairs = {'frag_0_frag_3': -2.5, 'frag_1_frag_3': -1.8, ...}
        >>> 
        >>> # Apply fp-CovaLED redistribution
        >>> fpled_pairs = compute_fpcovaled_interactions(
        ...     led_int, lig_frags, rec_frags, inter_pairs
        ... )
        >>> 
        >>> # Now sum(fpled_pairs.values()) == E_int
    """
    # Calculate intra-molecular contributions directly from LED matrix
    intra_ligand = 0.0
    for i, fpled_i in enumerate(ligand_fpled_frags):
        for fpled_j in ligand_fpled_frags[i:]:  # Upper triangle including diagonal
            try:
                val = led_int.loc[fpled_i, fpled_j]
                if pd.notna(val):
                    intra_ligand += float(val)
            except (KeyError, ValueError):
                pass
    
    # Calculate intra-molecular contributions for receptor
    intra_receptor = 0.0
    for i, fpled_i in enumerate(receptor_fpled_frags):
        for fpled_j in receptor_fpled_frags[i:]:  # Upper triangle including diagonal
            try:
                val = led_int.loc[fpled_i, fpled_j]
                if pd.notna(val):
                    intra_receptor += float(val)
            except (KeyError, ValueError):
                pass
    
    intra_total = intra_ligand + intra_receptor
    
    # Calculate sum of absolute values for scaling
    total_abs_inter = sum(abs(val) for val in fragment_pair_values.values())
    
    # Redistribute intra-molecular contributions proportionally
    if total_abs_inter > 0 and abs(intra_total) > 1e-6:
        fpled_interactions = {}
        for frag_pair, inter_val in fragment_pair_values.items():
            # Scale factor based on absolute magnitude
            scale_factor = abs(inter_val) / total_abs_inter
            redistribution = intra_total * scale_factor
            fpled_val = inter_val + redistribution
            fpled_interactions[frag_pair] = fpled_val
        
        return fpled_interactions
    else:
        # No redistribution needed (no intra-molecular contributions or no inter-molecular terms)
        return fragment_pair_values.copy()


def compute_led_interaction_matrix(
    orca_super: Union[str, Path],
    orca_ligand: Union[str, Path],
    orca_receptor: Union[str, Path],
    ligand_fpled_frags: Optional[List[int]] = None,
    receptor_fpled_frags: Optional[List[int]] = None,
    ligand_atom_indices: Optional[List[int]] = None,
    component: str = 'TOTAL',
    method: str = "DLPNO-CCSD(T)",
    conversion_factor: float = 2625.5,
    use_cache: bool = True,
    force_regenerate: bool = False,
    assign_mixed_fragments: str = 'error',
    verbose: bool = False
) -> pd.DataFrame:
    """
    High-level function to compute LED interaction matrix from ORCA output files.
    
    This function:
    1. Auto-detects fragment assignments from ORCA files (if not provided)
    2. Checks if LED Excel files exist, runs LEDAW if needed (vendored)
    3. Loads LED matrices (with optional caching for subsystems)
    4. Applies CovaLED Step 7 (subsystem subtraction)
    5. Returns the interaction LED matrix
    
    LEDAW Integration:
    LEDAW is bundled as a vendored dependency - no external installation required!
    Just import and use. LEDAW credit belongs to its original developers.
    
    For memory efficiency when processing many molecules:
    - Subsystem matrices are cached by default (set use_cache=False to disable)
    - Supersystem matrices are loaded fresh each time (unique per molecule)
    
    Args:
        orca_super: Path to supersystem ORCA .out file
        orca_ligand: Path to ligand subsystem ORCA .out file
        orca_receptor: Path to receptor subsystem ORCA .out file
        ligand_fpled_frags: List of fragment indices (1-based) belonging to ligand.
                           If None, auto-detected from ORCA file.
        receptor_fpled_frags: List of fragment indices (1-based) belonging to receptor.
                             If None, auto-detected from ORCA file.
        ligand_atom_indices: List of atom indices (0-based) that belong to ligand.
                            Required only if auto-detecting fragments.
                            If None, assumes first N/2 atoms are ligand.
        component: LED component to extract ('TOTAL', 'Electrostat', 'Exchange', etc.)
        method: QM method used ("DLPNO-CCSD(T)", "DLPNO-CCSD", "HFLD")
        conversion_factor: Energy conversion factor (default: 2625.5 for Hartree to kJ/mol)
        use_cache: Whether to cache subsystem matrices (default: True)
        force_regenerate: If True, regenerate LED files even if they exist
        assign_mixed_fragments: How to handle fragments spanning both molecules:
                               'error' (default), 'majority', 'ligand', 'receptor'
        verbose: If True, enable detailed LEDAW execution logging for debugging
    
    Returns:
        DataFrame with LED interaction matrix (same indexing as supersystem)
        
    Example:
        >>> # Auto-detect fragments (specify which atoms are ligand)
        >>> led_int = compute_led_interaction_matrix(
        ...     "complex/complex.out",
        ...     "ligand/ligand.out", 
        ...     "receptor/receptor.out",
        ...     ligand_atom_indices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        ... )
        >>> 
        >>> # Or manually specify fragment assignments
        >>> led_int = compute_led_interaction_matrix(
        ...     "complex/complex.out",
        ...     "ligand/ligand.out", 
        ...     "receptor/receptor.out",
        ...     ligand_fpled_frags=[1, 2, 3],
        ...     receptor_fpled_frags=[4, 5]
        ... )
        >>> 
        >>> # Get total interaction energy
        >>> e_int = extract_interaction_energy(led_int)
        >>> 
        >>> # Extract different components
        >>> elec_int = compute_led_interaction_matrix(
        ...     "complex/complex.out",
        ...     "ligand/ligand.out",
        ...     "receptor/receptor.out",
        ...     ligand_atom_indices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
        ...     component='Electrostat'
        ... )
    """
    # Auto-detect fragment assignments if not provided
    if ligand_fpled_frags is None or receptor_fpled_frags is None:
        atom_to_fragment, _ = parse_orca_fragments(orca_super)
        ligand_fpled_frags, receptor_fpled_frags = determine_fragment_groups(
            atom_to_fragment, ligand_atom_indices, assign_mixed=assign_mixed_fragments
        )
        print(f"  Auto-detected fragments: ligand={ligand_fpled_frags}, receptor={receptor_fpled_frags}")
    
    # Ensure LED Excel files exist (run LEDAW if needed)
    excel_super = _ensure_led_file(orca_super, method, conversion_factor, force_regenerate, verbose=verbose)
    excel_ligand = _ensure_led_file(orca_ligand, method, conversion_factor, force_regenerate, verbose=verbose)
    excel_receptor = _ensure_led_file(orca_receptor, method, conversion_factor, force_regenerate, verbose=verbose)
    
    # Load supersystem (unique per molecule, must load fresh)
    df_super = pd.read_excel(excel_super, sheet_name=component, index_col=0)
    
    # Load subsystems (optionally cached for efficiency)
    if use_cache:
        df_ligand = _get_cached_subsystem_led(excel_ligand, component)
        df_receptor = _get_cached_subsystem_led(excel_receptor, component)
    else:
        df_ligand = pd.read_excel(excel_ligand, sheet_name=component, index_col=0)
        df_receptor = pd.read_excel(excel_receptor, sheet_name=component, index_col=0)
    
    # Compute interaction matrix using CovaLED Step 7
    led_int = _compute_covaled_interaction_matrix(
        df_super, df_ligand, df_receptor,
        ligand_fpled_frags, receptor_fpled_frags
    )
    
    return led_int


def extract_interaction_energy(
    led_int: pd.DataFrame
) -> float:
    """
    Calculate total interaction energy from LED interaction matrix.
    
    Sums all elements in the LED interaction matrix to get E_int.
    After CovaLED Step 7, this equals the total interaction energy.
    
    Args:
        led_int: LED interaction matrix from compute_led_interaction_matrix
    
    Returns:
        Total interaction energy (in the same units as the LED matrix, typically kJ/mol)
        
    Example:
        >>> led_int = compute_led_interaction_matrix(...)
        >>> e_int = extract_interaction_energy(led_int)
        >>> print(f"Interaction energy: {e_int:.4f} kJ/mol")
    """
    return float(led_int.sum().sum())


def _process_fp_led_matrices_in_memory(
    standard_led_matrices: Dict[str, pd.DataFrame],
    method: str,
    main_filenames: List[Union[str, Path]],
    alternative_filenames: List[Union[str, Path]],
    conversion_factor: float
) -> Dict[str, pd.DataFrame]:
    """
    Process fp-LED matrices in memory from standard LED matrices.
    
    This is an in-memory version of process_fp_LED_matrices that works with DataFrames
    instead of Excel files.
    """
    from xbpy._vendor.ledaw_package.nbody_engine import (
        compute_fp_el_prep,
        compute_bulk_solvation_contribution
    )
    
    # Process REF sheet
    df_ref = standard_led_matrices.get('REF', pd.DataFrame())
    if df_ref.empty:
        raise ValueError("REF matrix not found in standard LED matrices")
    
    df_ref.columns = [int(c) for c in df_ref.columns]
    df_ref.index = [int(i) for i in df_ref.index]
    
    df_ref_distributed = compute_fp_el_prep(df_ref)
    df_ref_distributed.columns.name = 'REF-EL-PREP'
    
    # Get SOLV values and calculate SOLV interaction energy
    solv_values, e_solv_int_en, _, _, _ = compute_bulk_solvation_contribution(
        main_filenames, alternative_filenames, conversion_factor, method
    )
    
    fp_led_matrices = {}
    
    # Load Electrostat and Exchange sheets for all methods
    df_electrostat = standard_led_matrices.get('Electrostat', pd.DataFrame())
    df_exchange = standard_led_matrices.get('Exchange', pd.DataFrame())
    df_electrostat.columns = [int(c) for c in df_electrostat.columns]
    df_electrostat.index = [int(i) for i in df_electrostat.index]
    df_exchange.columns = [int(c) for c in df_exchange.columns]
    df_exchange.index = [int(i) for i in df_exchange.index]
    
    if method.lower() in ["dlpno-ccsd(t)", "dlpno-ccsd"]:
        # Read TOTAL
        df_total = standard_led_matrices.get('TOTAL', pd.DataFrame())
        df_total.columns = [int(c) for c in df_total.columns]
        df_total.index = [int(i) for i in df_total.index]
        
        # Subtract SOLV to get standard LED total without SOLV
        df_solv = standard_led_matrices.get('SOLV', pd.DataFrame())
        if not df_solv.empty:
            df_solv.columns = [int(c) for c in df_solv.columns]
            df_solv.index = [int(i) for i in df_solv.index]
            std_df_total_no_solv = df_total - df_solv
        else:
            std_df_total_no_solv = df_total
        
        # Correlation EL-PREP on standard total without SOLV
        df_total_distributed = compute_fp_el_prep(std_df_total_no_solv)
        df_total_distributed.columns.name = 'C-CCSD(T)-EL-PREP' if method.lower() == "dlpno-ccsd(t)" else 'C-CCSD-EL-PREP'
        
        df_c_ccsd = df_total_distributed - df_ref_distributed
        df_c_ccsd.columns.name = 'CCSD(T)-EL-PREP' if method.lower() == "dlpno-ccsd(t)" else 'CCSD-EL-PREP'
        
        # DISP - INTER-NODISP decomposition
        disp_sheet = 'Disp CCSD' if method.lower() == "dlpno-ccsd" else 'Disp CCSD(T)'
        inter_nondisp_sheet = 'Inter-NonDisp-C-CCSD' if method.lower() == "dlpno-ccsd" else 'Inter-NonDisp-C-CCSD(T)'
        
        df_disp_ccsd = standard_led_matrices.get(disp_sheet, pd.DataFrame())
        df_inter_nondisp = standard_led_matrices.get(inter_nondisp_sheet, pd.DataFrame())
        df_disp_ccsd.columns = [int(c) for c in df_disp_ccsd.columns]
        df_disp_ccsd.index = [int(i) for i in df_disp_ccsd.index]
        df_inter_nondisp.columns = [int(c) for c in df_inter_nondisp.columns]
        df_inter_nondisp.index = [int(i) for i in df_inter_nondisp.index]
        
        # Reconstruct total without SOLV from fp-LED logic
        df_ref_final = df_electrostat + df_exchange + df_ref_distributed
        df_c_ccsd_final = df_c_ccsd + df_disp_ccsd + df_inter_nondisp
        df_total_no_solv = df_ref_final + df_c_ccsd_final
        
        # Apply fp-LED SOLV if needed
        if e_solv_int_en != 0:
            solv_matrix = df_total_no_solv * (e_solv_int_en / df_total_no_solv.sum().sum())
            df_total_final = df_total_no_solv + solv_matrix
            fp_led_matrices['SOLV'] = solv_matrix
        else:
            df_total_final = df_total_no_solv
        
        fp_led_matrices['TOTAL'] = df_total_final
        fp_led_matrices['REF'] = df_ref_final
        fp_led_matrices['Electrostat'] = df_electrostat
        fp_led_matrices['Exchange'] = df_exchange
        fp_led_matrices['REF-EL-PREP'] = df_ref_distributed
        fp_led_matrices['C-CCSD(T)' if method.lower() == "dlpno-ccsd(t)" else 'C-CCSD'] = df_c_ccsd_final
        fp_led_matrices[disp_sheet] = df_disp_ccsd
        fp_led_matrices[inter_nondisp_sheet] = df_inter_nondisp
        fp_led_matrices['C-CCSD(T)-EL-PREP' if method.lower() == "dlpno-ccsd(t)" else 'C-CCSD-EL-PREP'] = df_c_ccsd
        fp_led_matrices['CCSD(T)-EL-PREP' if method.lower() == "dlpno-ccsd(t)" else 'CCSD-EL-PREP'] = df_total_distributed
        
    elif method.lower() == "hfld":
        # Calculate REF = Electrostat + Exchange + REF-EL-PREP
        df_ref_final = df_electrostat + df_exchange + df_ref_distributed
        
        # Calculate TOTAL = REF + Disp HFLD
        df_disp_hfld = standard_led_matrices.get('Disp HFLD', pd.DataFrame())
        df_disp_hfld.columns = [int(c) for c in df_disp_hfld.columns]
        df_disp_hfld.index = [int(i) for i in df_disp_hfld.index]
        
        df_total_hfld = df_ref_final + df_disp_hfld
        df_total_hfld.columns.name = 'TOTAL'
        
        # If SOLV interaction energy is not zero, calculate and add the SOLV matrix
        if e_solv_int_en != 0:
            solv_matrix = df_total_hfld * (e_solv_int_en / df_total_hfld.sum().sum())
            df_total_hfld += solv_matrix
            fp_led_matrices['SOLV'] = solv_matrix
        
        fp_led_matrices['TOTAL'] = df_total_hfld
        fp_led_matrices['REF'] = df_ref_final
        fp_led_matrices['Electrostat'] = df_electrostat
        fp_led_matrices['Exchange'] = df_exchange
        fp_led_matrices['REF-EL-PREP'] = df_ref_distributed
        fp_led_matrices['Disp HFLD'] = df_disp_hfld
    
    return fp_led_matrices


def compute_fp_led_interactions_with_mapping(
    main_filenames: List[Union[str, Path]],
    alternative_filenames: List[Union[str, Path]],
    fragment_mappings: Dict[str, Dict[int, int]],
    conversion_factor: float = 2625.5,
    method: str = "DLPNO-CCSD(T)",
    LEDAW_output_path: Optional[Union[str, Path]] = None,
    use_temp_dir: bool = True,
    use_dataframes: bool = False,
    cleanup_excel: bool = True,
    verbose: bool = False
) -> Dict[str, float]:
    """
    Compute fp-LED pairwise interactions using custom fragment mappings.
    
    This function allows you to specify fragment mappings directly instead of
    matching by position. Returns fp-LED pairwise interactions as a dictionary.
    
    Args:
        main_filenames: List of ORCA output file paths [supersystem, subsystem1, subsystem2, ...]
        alternative_filenames: List of alternative ORCA output file paths (same length as main_filenames)
        fragment_mappings: Dictionary mapping system labels to fragment mappings.
                         Format: {
                             "SUPERSYS": {frag_id: frag_id, ...},  # Identity mapping for supersystem
                             "SUBSYS1": {subsys_frag_id: supersys_frag_id, ...},
                             "SUBSYS2": {subsys_frag_id: supersys_frag_id, ...},
                             ...
                         }
                         Fragment IDs are 1-based (as in ORCA output).
        conversion_factor: Energy conversion factor (default: 2625.5 for Hartree to kJ/mol)
        method: QM method used ("DLPNO-CCSD(T)", "DLPNO-CCSD", "HFLD")
        LEDAW_output_path: Directory for LEDAW output (default: same as first file's directory)
        use_temp_dir: If True, use a temporary directory for Excel files and clean up after (default: True)
        use_dataframes: If True, process fp-LED matrices in memory using DataFrames instead of Excel (default: False)
        cleanup_excel: If True, delete Excel files after extracting data (default: True, only used if use_temp_dir=False)
        verbose: If True, print detailed progress information
    
    Returns:
        Dictionary of fp-LED pairwise interactions in format:
        {"frag_i_frag_j": energy, ...} where i <= j (upper triangle only)
        Energies are in the same units as conversion_factor (typically kJ/mol)
        
    Example:
        >>> mappings = {
        ...     "SUPERSYS": {1: 1, 2: 2, 3: 3, 4: 4},
        ...     "SUBSYS1": {1: 1, 2: 2, 3: 3},  # Fragments 1,2,3 map to supersystem 1,2,3
        ...     "SUBSYS2": {1: 4}  # Fragment 1 maps to supersystem 4
        ... }
        >>> fp_interactions = compute_fp_led_interactions_with_mapping(
        ...     ["complex.out", "ligand.out", "receptor.out"],
        ...     ["", "", ""],
        ...     mappings,
        ...     use_temp_dir=True  # Use temporary directory, automatically cleaned up
        ... )
        >>> # fp_interactions = {"1_1": -5.234, "1_2": -2.156, "1_4": -1.234, ...}
    """
    import tempfile
    import shutil
    from xbpy._vendor.ledaw_package.nbody_engine import (
        construct_label_mappings as _original_construct_label_mappings,
        normalize_path
    )
    
    if not _LEDAW_AVAILABLE:
        raise ImportError(f"LEDAW not available: {_LEDAW_IMPORT_ERROR}")
    
    # Set up output path - use temporary directory if requested
    temp_dir = None
    if use_temp_dir:
        temp_dir = tempfile.mkdtemp(prefix="ledaw_")
        LEDAW_output_path = temp_dir
        if verbose:
            print(f"Using temporary directory for LEDAW output: {LEDAW_output_path}")
    else:
        if LEDAW_output_path is None:
            LEDAW_output_path = Path(main_filenames[0]).parent
        else:
            LEDAW_output_path = Path(LEDAW_output_path)
    
    LEDAW_output_path = normalize_path(str(LEDAW_output_path))
    
    # Create a custom construct_label_mappings function that uses provided mappings
    def _custom_construct_label_mappings(main_fns, alt_fns, output_path):
        """Custom version that uses provided fragment mappings."""
        # Convert fragment mappings to the format expected by LEDAW
        main_label_mappings_ghost_free = {}
        alt_label_mappings_ghost_free = {}
        
        # SUPERSYS mapping (should be identity)
        if "SUPERSYS" in fragment_mappings:
            main_label_mappings_ghost_free["SUPERSYS"] = fragment_mappings["SUPERSYS"].copy()
        else:
            # Default: identity mapping
            from xbpy._vendor.ledaw_package.nbody_engine import extract_fragments
            supersystem_file = normalize_path(main_fns[0])
            supersystem_frags, super_ghost_flags, _ = extract_fragments(supersystem_file)
            main_label_mappings_ghost_free["SUPERSYS"] = {
                k: k for k, v in super_ghost_flags.items() if not v
            }
        
        alt_label_mappings_ghost_free["SUPERSYS"] = main_label_mappings_ghost_free["SUPERSYS"].copy()
        
        # SUBSYS mappings
        for i in range(1, len(main_fns)):
            subsys_label = f"SUBSYS{i}"
            if subsys_label in fragment_mappings:
                main_label_mappings_ghost_free[subsys_label] = fragment_mappings[subsys_label].copy()
            else:
                # Default: empty mapping
                main_label_mappings_ghost_free[subsys_label] = {}
            
            alt_label_mappings_ghost_free[subsys_label] = main_label_mappings_ghost_free[subsys_label].copy()
        
        # Determine BSSE (check if any subsystem has ghost atoms)
        bsse_found = False
        for i in range(1, len(main_fns)):
            from xbpy._vendor.ledaw_package.nbody_engine import extract_fragments
            subsystem_file = normalize_path(main_fns[i])
            _, ghost_flags, bsse = extract_fragments(subsystem_file)
            bsse_found |= bsse
        
        if verbose:
            print("======== Fragment Mappings (Custom) ========")
            print("main_label_mappings_ghost_free:", main_label_mappings_ghost_free)
            print("alt_label_mappings_ghost_free:", alt_label_mappings_ghost_free)
            print("bsse_found:", bsse_found)
        
        return main_label_mappings_ghost_free, alt_label_mappings_ghost_free, bsse_found
    
    # Temporarily replace construct_label_mappings
    import xbpy._vendor.ledaw_package.nbody_engine as nbody_engine
    original_func = nbody_engine.construct_label_mappings
    nbody_engine.construct_label_mappings = _custom_construct_label_mappings
    
    try:
        # Run LEDAW engine
        if verbose:
            print(f"Running LEDAW with custom fragment mappings...")
        
        _LEDAW_ENGINE(
            main_filenames=[str(f) for f in main_filenames],
            alternative_filenames=[str(f) if f else "" for f in alternative_filenames],
            conversion_factor=conversion_factor,
            method=method,
            LEDAW_output_path=LEDAW_output_path,
            relabel_mapping=None,
            use_ref_as_rhf_in_hfld=None
        )
        
        # Read standard LED matrices from Excel (needed for fp-LED processing)
        standard_excel_path = Path(LEDAW_output_path) / "Summary_Standard_LED_matrices.xlsx"
        if not standard_excel_path.exists():
            raise FileNotFoundError(f"Standard LED Excel file not found: {standard_excel_path}")
        
        # Load all standard LED matrices into memory
        standard_led_matrices = {}
        with pd.ExcelFile(standard_excel_path) as xl:
            for sheet_name in xl.sheet_names:
                df = pd.read_excel(xl, sheet_name=sheet_name, index_col=0)
                standard_led_matrices[sheet_name] = df
        
        # Process fp-LED matrices in memory (no Excel files for fp-LED)
        if use_dataframes:
            # Process entirely in memory
            fp_led_matrices = _process_fp_led_matrices_in_memory(
                standard_led_matrices, method, main_filenames, alternative_filenames, conversion_factor
            )
            df_fp_led = fp_led_matrices.get('TOTAL', pd.DataFrame())
        else:
            # Use existing Excel-based approach (backward compatibility)
            fp_excel_path = Path(LEDAW_output_path) / "Summary_fp-LED_matrices.xlsx"
            if not fp_excel_path.exists():
                raise FileNotFoundError(f"fp-LED Excel file not found: {fp_excel_path}")
            df_fp_led = pd.read_excel(fp_excel_path, sheet_name='TOTAL', index_col=0)
        
        # Convert to dictionary format
        fp_interactions = {}
        for frag_i in df_fp_led.index:
            if pd.notna(frag_i):
                for frag_j in df_fp_led.columns:
                    if pd.notna(frag_j):
                        # Get the interaction energy for this pair
                        # Use upper triangle to avoid duplicates (LED matrices are symmetric)
                        if frag_i <= frag_j:
                            energy = df_fp_led.loc[frag_i, frag_j]
                            if pd.notna(energy):
                                key = f"{int(frag_i)}_{int(frag_j)}"
                                fp_interactions[key] = float(energy)
        
        # Clean up Excel files if requested
        if use_temp_dir:
            # If using temp directory, always clean it up
            try:
                shutil.rmtree(temp_dir)
                if verbose:
                    print(f"Cleaned up temporary directory: {temp_dir}")
            except Exception as e:
                if verbose:
                    print(f"Error cleaning up temporary directory: {e}")
        elif cleanup_excel:
            # Only clean up Excel files if not using temp dir
            try:
                import glob
                excel_files = glob.glob(str(Path(LEDAW_output_path) / "*.xlsx"))
                for excel_file in excel_files:
                    try:
                        os.remove(excel_file)
                        if verbose:
                            print(f"Deleted: {excel_file}")
                    except Exception as e:
                        if verbose:
                            print(f"Could not delete {excel_file}: {e}")
            except Exception as e:
                if verbose:
                    print(f"Error during cleanup: {e}")
        
        return fp_interactions
        
    finally:
        # Restore original function
        nbody_engine.construct_label_mappings = original_func
        # Ensure temp directory is cleaned up even on error
        if use_temp_dir and temp_dir and os.path.exists(temp_dir):
            try:
                shutil.rmtree(temp_dir)
            except Exception:
                pass


def compute_fp_led_interactions(
    orca_out_path: Union[str, Path],
    method: str = "DLPNO-CCSD(T)",
    conversion_factor: float = 2625.5,
    force_regenerate: bool = False,
    verbose: bool = False
) -> Dict[str, float]:
    """
    Compute fp-LED pairwise interactions for a single ORCA output file.
    
    This is a convenience function that runs LEDAW and extracts fp-LED interactions
    for a single system (no subsystem subtraction).
    
    Args:
        orca_out_path: Path to ORCA .out file
        method: QM method used ("DLPNO-CCSD(T)", "DLPNO-CCSD", "HFLD")
        conversion_factor: Energy conversion factor (default: 2625.5 for Hartree to kJ/mol)
        force_regenerate: If True, regenerate LED files even if they exist
        verbose: If True, enable detailed output
    
    Returns:
        Dictionary of fp-LED pairwise interactions in format:
        {"frag_i_frag_j": energy, ...} where i <= j (upper triangle only)
        
    Example:
        >>> fp_interactions = compute_fp_led_interactions("molecule.out")
        >>> # fp_interactions = {"1_1": -5.234, "1_2": -2.156, "2_2": -3.421, ...}
    """
    orca_out_path = Path(orca_out_path)
    
    # Ensure LED Excel file exists (run LEDAW if needed)
    excel_path = _ensure_led_file(
        orca_out_path, 
        method=method, 
        conversion_factor=conversion_factor,
        force_regenerate=force_regenerate,
        verbose=verbose
    )
    
    # Read fp-LED results
    fp_excel_path = excel_path.parent / "Summary_fp-LED_matrices.xlsx"
    if not fp_excel_path.exists():
        raise FileNotFoundError(f"fp-LED Excel file not found: {fp_excel_path}. LEDAW may not have generated fp-LED matrices.")
    
    # Load fp-LED TOTAL matrix
    df_fp_led = pd.read_excel(fp_excel_path, sheet_name='TOTAL', index_col=0)
    
    # Convert to dictionary format
    fp_interactions = {}
    for frag_i in df_fp_led.index:
        if pd.notna(frag_i):
            for frag_j in df_fp_led.columns:
                if pd.notna(frag_j):
                    # Get the interaction energy for this pair
                    # Use upper triangle to avoid duplicates (LED matrices are symmetric)
                    if frag_i <= frag_j:
                        energy = df_fp_led.loc[frag_i, frag_j]
                        if pd.notna(energy):
                            key = f"{int(frag_i)}_{int(frag_j)}"
                            fp_interactions[key] = float(energy)
    
    return fp_interactions

