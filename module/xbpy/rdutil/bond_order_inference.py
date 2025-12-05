
from ..mathutils.geometry import calculate_angle
from itertools import combinations
from .geometry import position
from .util import possible_geometries, ideal_bond_lengths
import numpy as np
from rdkit import Chem

def correct_bond_orders(mol, add_hydrogens = False, try_full = True):
    try:
        from z3 import Optimize, sat, Solver, IntVal
    except ImportError:
        raise ImportError("The z3 package is required for this function. Try to install via 'sudo apt install z3' and 'pip install z3-solver'.")

    all_bond_length = {}
    for bond in mol.GetBonds():
        length = np.linalg.norm(position(bond.GetBeginAtom()) - position(bond.GetEndAtom()))
        all_bond_length[(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())] = length
        all_bond_length[(bond.GetEndAtomIdx(), bond.GetBeginAtomIdx())] = length
    all_pairwise_bond_angles = _get_all_pairwise_bond_angles(mol)
    violating_atoms = _get_violating_atom_clusters(mol, all_pairwise_bond_angles, add_hydrogens = add_hydrogens)
    if len(violating_atoms) == 0:
        return mol
    # Expand the set by one-hop neighbors so substructures include bonded context (incl. hydrogens)
    expanded_atoms = set(violating_atoms)
    for idx in list(violating_atoms):
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            expanded_atoms.add(nbr.GetIdx())
    connected_subgraphs = _get_connected_subgraphs(mol, expanded_atoms)
    new_mol = Chem.RWMol(mol)
    try:
        for subgraph in connected_subgraphs:
            # Print substructure atoms concisely (debug only)
            import os
            if os.environ.get('XB_DEBUG_SAT') == '1':
                print("Substructure:", [f"{mol.GetAtomWithIdx(i).GetSymbol()}{i}" for i in sorted(subgraph)])
            
            # First try to soften constraints iteratively until SAT
            softened_atoms, strict_atoms = _resolve_unsat_cores_iteratively(mol, subgraph, all_pairwise_bond_angles, all_bond_length)
            
            if not strict_atoms:
                pass
            
            # Now solve with all atoms, but with softened constraints for problematic ones
            
            formula, objective, variable_dict, named_constraints = _construct_SMT_formulation(mol, subgraph, all_pairwise_bond_angles, all_bond_length, excluded_atoms=softened_atoms)
            
            # First check if the system is SAT at all
            s = Solver()
            # Track all named constraints only if debugging
            import os
            if os.environ.get('XB_DEBUG_SAT') == '1':
                for cname, cexpr in named_constraints.items():
                    s.assert_and_track(cexpr, cname)
            sat_result = s.check()
            
            if sat_result != sat:
                raise ValueError(f"System is not satisfiable even after constraint softening. SAT result: {sat_result}")
            
            # Debug: dump SAT model (pre-optimization) gated by env var XB_DEBUG_SAT=1
            import os
            if os.environ.get('XB_DEBUG_SAT') == '1':
                try:
                    sat_model = s.model()
                    print("\n--- DEBUG: Pre-optimization SAT model (non-zero selections) ---")
                    # Bond orders selected (==1)
                    bond_selected = []
                    charge_selected = []
                    geom_selected = []
                    for key, var in variable_dict.items():
                        try:
                            val = sat_model.evaluate(var, model_completion=True)
                            if hasattr(val, 'as_long'):
                                ival = val.as_long()
                            else:
                                ival = int(str(val))
                        except Exception:
                            continue
                        if ival != 1:
                            continue
                        if isinstance(key, tuple) and len(key) == 3:
                            a1, a2, bo = key
                            bond_selected.append((a1, a2, bo))
                        elif isinstance(key, tuple) and len(key) == 2:
                            a, ch = key
                            charge_selected.append((a, ch))
                        elif isinstance(key, tuple) and len(key) == 4:
                            a, nc, ch, ang = key
                            geom_selected.append((a, nc, ch, ang))
                    if bond_selected:
                        print("Bond orders (==1):")
                        for a1, a2, bo in bond_selected:
                            print(f"  {mol.GetAtomWithIdx(a1).GetSymbol()}{a1}-{mol.GetAtomWithIdx(a2).GetSymbol()}{a2}: {bo}")
                    if charge_selected:
                        print("Formal charges (==1):")
                        for a, ch in charge_selected:
                            ctype = "softened" if a in softened_atoms else "strict"
                            print(f"  {mol.GetAtomWithIdx(a).GetSymbol()}{a}: {ch} ({ctype})")
                    if geom_selected:
                        print("Geometries (==1):")
                        for a, nc, ch, ang in geom_selected:
                            print(f"  {mol.GetAtomWithIdx(a).GetSymbol()}{a}: neighbors={nc}, charge={ch}, angles={ang}")
                except Exception as e:
                    print(f"(debug model dump failed: {e})")

            # If SAT, use Optimize to find the best solution
            
            opt = Optimize()
            opt.add(formula)
            opt.minimize(objective)
            opt_result = opt.check()
            
            if opt_result == sat:
                model = opt.model()
                new_bond_orders, assigned_charges = _decode_solution(model)
                
                # Print only assignments relevant to softened atoms (only if any softened)
                if softened_atoms:
                    print(f"\n=== FINAL ASSIGNMENTS (softened atoms only) ===")
                    # Bonds per softened atom
                    any_bond_printed = False
                    for a in sorted(softened_atoms):
                        atom_symbol = mol.GetAtomWithIdx(a).GetSymbol()
                        printed_for_atom = False
                        for (a1, a2), bo in new_bond_orders.items():
                            if a1 == a or a2 == a:
                                if not any_bond_printed:
                                    print("Bond orders:")
                                    any_bond_printed = True
                                nbr = a2 if a1 == a else a1
                                nbr_symbol = mol.GetAtomWithIdx(nbr).GetSymbol()
                                print(f"  {atom_symbol}{a}-{nbr_symbol}{nbr}: {bo}")
                                printed_for_atom = True
                        if not printed_for_atom:
                            print(f"  {atom_symbol}{a}: no bond order selection (all 0)")
                    # Charges for softened atoms
                    any_charge_printed = False
                    for a in sorted(softened_atoms):
                        atom_symbol = mol.GetAtomWithIdx(a).GetSymbol()
                        if a in assigned_charges:
                            if not any_charge_printed:
                                print("Formal charges:")
                                any_charge_printed = True
                            print(f"  {atom_symbol}{a}: {assigned_charges[a]} (softened)")
                    if not any_charge_printed:
                        print("Formal charges: none selected for softened atoms")
                
                # Apply the assignments to the molecule
                for (atom1, atom2), bond_order in new_bond_orders.items():
                    if bond_order == 0.0:
                        new_mol.RemoveBond(int(atom1), int(atom2))
                    else:
                        new_mol.GetBondBetweenAtoms(int(atom1), int(atom2)).SetBondType(Chem.BondType(bond_order))
                for atom, charge in assigned_charges.items():
                    new_mol.GetAtomWithIdx(atom).SetFormalCharge(charge)
                
            else:
                raise ValueError(f"Optimization failed for atoms {[mol.GetAtomWithIdx(idx).GetSymbol() + ':' + str(idx) for idx in subgraph]}. Result: {opt_result}")
        return new_mol
    except ValueError as e:
        if try_full:
            
            # Try to soften constraints iteratively for the full molecule
            softened_atoms, strict_atoms = _resolve_unsat_cores_iteratively(mol, range(mol.GetNumAtoms()), all_pairwise_bond_angles, all_bond_length)
            
            if not strict_atoms:
                pass
            
            # Now solve with all atoms, but with softened constraints for problematic ones
            
            formula, objective, variable_dict, named_constraints = _construct_SMT_formulation(mol, range(mol.GetNumAtoms()), all_pairwise_bond_angles, all_bond_length, excluded_atoms=softened_atoms)
            
            # First check if the system is SAT at all
            s = Solver()
            s.add(formula)
            sat_result = s.check()
            
            if sat_result != sat:
                raise ValueError(f"System is not satisfiable even after constraint softening. SAT result: {sat_result}")
            
            # Debug: dump SAT model (pre-optimization) gated by env var XB_DEBUG_SAT=1
            import os
            if os.environ.get('XB_DEBUG_SAT') == '1':
                try:
                    sat_model = s.model()
                    print("\n--- DEBUG: Pre-optimization SAT model (non-zero selections) ---")
                    bond_selected = []
                    charge_selected = []
                    geom_selected = []
                    for key, var in variable_dict.items():
                        try:
                            val = sat_model.evaluate(var, model_completion=True)
                            if hasattr(val, 'as_long'):
                                ival = val.as_long()
                            else:
                                ival = int(str(val))
                        except Exception:
                            continue
                        if ival != 1:
                            continue
                        if isinstance(key, tuple) and len(key) == 3:
                            a1, a2, bo = key
                            bond_selected.append((a1, a2, bo))
                        elif isinstance(key, tuple) and len(key) == 2:
                            a, ch = key
                            charge_selected.append((a, ch))
                        elif isinstance(key, tuple) and len(key) == 4:
                            a, nc, ch, ang = key
                            geom_selected.append((a, nc, ch, ang))
                    if bond_selected:
                        print("Bond orders (==1):")
                        for a1, a2, bo in bond_selected:
                            print(f"  {mol.GetAtomWithIdx(a1).GetSymbol()}{a1}-{mol.GetAtomWithIdx(a2).GetSymbol()}{a2}: {bo}")
                    if charge_selected:
                        print("Formal charges (==1):")
                        for a, ch in charge_selected:
                            ctype = "softened" if a in softened_atoms else "strict"
                            print(f"  {mol.GetAtomWithIdx(a).GetSymbol()}{a}: {ch} ({ctype})")
                    if geom_selected:
                        print("Geometries (==1):")
                        for a, nc, ch, ang in geom_selected:
                            print(f"  {mol.GetAtomWithIdx(a).GetSymbol()}{a}: neighbors={nc}, charge={ch}, angles={ang}")
                except Exception as e:
                    print(f"(debug model dump failed: {e})")

            # If SAT, use Optimize to find the best solution
            
            opt = Optimize()
            opt.add(formula)
            opt.minimize(objective)
            opt_result = opt.check()
            
            if opt_result == sat:
                model = opt.model()
                new_bond_orders, assigned_charges = _decode_solution(model)
                
                # Print only assignments relevant to softened atoms (only if any softened)
                if softened_atoms:
                    print(f"\n=== FINAL ASSIGNMENTS (softened atoms only) ===")
                    soft_bonds = [
                        (a1, a2, bo) for (a1, a2), bo in new_bond_orders.items()
                        if a1 in softened_atoms or a2 in softened_atoms
                    ]
                    if soft_bonds:
                        print(f"Bond orders:")
                        for a1, a2, bo in soft_bonds:
                            atom1_symbol = mol.GetAtomWithIdx(a1).GetSymbol()
                            atom2_symbol = mol.GetAtomWithIdx(a2).GetSymbol()
                            print(f"  {atom1_symbol}{a1}-{atom2_symbol}{a2}: {bo}")
                    soft_charges = [(a, ch) for a, ch in assigned_charges.items() if a in softened_atoms]
                    if soft_charges:
                        print(f"Formal charges:")
                        for atom, charge in soft_charges:
                            atom_symbol = mol.GetAtomWithIdx(atom).GetSymbol()
                            print(f"  {atom_symbol}{atom}: {charge} (softened)")
                
                # Apply the assignments to the molecule
                for (atom1, atom2), bond_order in new_bond_orders.items():
                    if bond_order == 0.0:
                        new_mol.RemoveBond(int(atom1), int(atom2))
                    else:
                        new_mol.GetBondBetweenAtoms(int(atom1), int(atom2)).SetBondType(Chem.BondType(bond_order))
                for atom, charge in assigned_charges.items():
                    new_mol.GetAtomWithIdx(atom).SetFormalCharge(charge)
                return new_mol
            else:
                raise ValueError(f"Optimization failed for the full molecule even after removing problematic atoms. Result: {opt_result}")


def _decode_solution(model):
    new_bond_orders = {}
    assigned_charges = {}
    for decl in model.decls():
        # Get the variable name
        key = decl.name()
        # Retrieve the value assigned to the variable
        value = model[decl]
        if value.as_long() == 1:
            if key.startswith("bond_order"):
                _, _, atom1, atom2, bond_order = key.split("_")
                new_bond_orders[(int(atom1), int(atom2))] = float(bond_order)
            elif key.startswith("charge"):
                _, atom, charge = key.split("_")
                assigned_charges[int(atom)] = int(charge)
    return new_bond_orders, assigned_charges

# .1 A deviation is as bad as 10 degrees
angle_deviation_weighting = 1
# might want to expand considered atoms by one layer of neighbors

ideal_bond_length_order_indices = {1.0: 0, 2.0: 2, 3.0: 3}

def _construct_SMT_formulation(mol, indices, all_pairwise_bond_angles, all_bond_Lengths, angle_tolerance=20, enable_atom_activation=False, excluded_atoms=None):
    try:
        from z3 import Or, And, Int, RealVal, Sum, Solver, Implies, IntVal, Bool, BoolVal
    except ImportError:
        raise ImportError("The z3 package is required for this function. Try to install via 'sudo apt install z3' and 'pip install z3-solver'.")

    considered_atoms = [mol.GetAtomWithIdx(idx) for idx in indices]
    penalty_expressions = [RealVal(0)]
    
    # Atom activation variables for unsat core analysis
    atom_activation_vars = {}
    if enable_atom_activation:
        for atom_idx in indices:
            atom_activation_vars[atom_idx] = Bool(f"atom_active_{atom_idx}")
    
    # Set of excluded atoms (will have relaxed constraints)
    excluded_atoms = excluded_atoms or set()
    
    # Track constraints with meaningful names based on atoms
    named_constraints = {}
    
    # variables for each bond order
    bond_order_variables = {}
    bond_order_exclusivity = []
    bond_order_constraints = []
    bond_order_is_bool = []

    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in indices and bond.GetEndAtomIdx() in indices:
            # bond order = 0 => no bond between the atoms
            for bond_order in range(4):
                bond_order = float(bond_order)
                symbol = Int(f"bond_order_{bond.GetBeginAtomIdx()}_{bond.GetEndAtomIdx()}_{bond_order}")
                bond_order_variables[(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond_order)] = symbol
                bond_order_variables[(bond.GetEndAtomIdx(), bond.GetBeginAtomIdx(), bond_order)] = symbol
                # Booleanity guarded/relaxed
                b_atom1 = bond.GetBeginAtomIdx()
                b_atom2 = bond.GetEndAtomIdx()
                if b_atom1 in excluded_atoms or b_atom2 in excluded_atoms:
                    bool_constraint = BoolVal(True)
                    bool_name = f"bool_bond_relaxed_{b_atom1}_{b_atom2}_{int(bond_order)}"
                else:
                    base_bool = Or(symbol == 0, symbol == 1)
                    if enable_atom_activation:
                        bool_constraint = Implies(And(atom_activation_vars[b_atom1], atom_activation_vars[b_atom2]), base_bool)
                        bool_name = f"bool_bond_{b_atom1}_{b_atom2}_{int(bond_order)}_conditional"
                    else:
                        bool_constraint = base_bool
                        bool_name = f"bool_bond_{b_atom1}_{b_atom2}_{int(bond_order)}"
                bond_order_is_bool.append(bool_constraint)
                named_constraints[bool_name] = bool_constraint
                # add bond length deviation penalty
                if bond_order > 0:
                    bond_length = all_bond_Lengths[(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())]
                    ideal_bond_length = ideal_bond_lengths[bond.GetBeginAtom().GetSymbol()][bond.GetEndAtom().GetSymbol()][ideal_bond_length_order_indices[bond_order]]
                    penalty_expressions.append((ideal_bond_length - bond_length) * symbol)
                else:
                    penalty_expressions.append(10 * symbol)

            # bond order exclusivity
            atom1_idx = bond.GetBeginAtomIdx()
            atom2_idx = bond.GetEndAtomIdx()
            
            # Check if either atom is excluded
            atom1_excluded = atom1_idx in excluded_atoms
            atom2_excluded = atom2_idx in excluded_atoms
            
            if atom1_excluded or atom2_excluded:
                # For excluded atoms, just require non-negative values (relaxed constraint)
                relaxed_constraint = And([bond_order_variables[(atom1_idx, atom2_idx, bond_order)] >= 0 for bond_order in range(4)])
                constraint = relaxed_constraint
                constraint_name = f"bond_relaxed_{atom1_idx}_{atom2_idx}"
            else:
                # For non-excluded atoms, enforce strict exclusivity
                base_constraint = Sum([bond_order_variables[(atom1_idx, atom2_idx, bond_order)] for bond_order in range(4)]) == 1
                
                if enable_atom_activation:
                    # Only enforce this constraint if both atoms are active
                    constraint = Implies(And(atom_activation_vars[atom1_idx], atom_activation_vars[atom2_idx]), base_constraint)
                    constraint_name = f"bond_exclusivity_{atom1_idx}_{atom2_idx}_conditional"
                else:
                    constraint = base_constraint
                    constraint_name = f"bond_exclusivity_{atom1_idx}_{atom2_idx}"
            
            bond_order_exclusivity.append(constraint)
            named_constraints[constraint_name] = constraint
        elif bond.GetBeginAtomIdx() in indices or bond.GetEndAtomIdx() in indices:
            for bond_order in range(4):
                bond_order = float(bond_order)
                bond_order_variables[(bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx(), bond_order)] = IntVal(int(bond.GetBondTypeAsDouble() == bond_order))
                bond_order_variables[(bond.GetEndAtom().GetIdx(), bond.GetBeginAtom().GetIdx(), bond_order)] = IntVal(int(bond.GetBondTypeAsDouble() == bond_order))

    # determine possible charges for each element
    possible_charges = {}
    for atom in considered_atoms:
        this_possible_charges = possible_charges.get(atom.GetSymbol(), set())
        for (neighbor_count, atom_charge) in possible_geometries[atom.GetSymbol()].keys():
            this_possible_charges.add(atom_charge)
        possible_charges[atom.GetSymbol()] = this_possible_charges

    # variable for each charge
    charge_variables = {}
    charge_exclusivity = []
    charge_is_bool = []
    for atom in considered_atoms:
        for possible_charge in possible_charges[atom.GetSymbol()]:
            symbol = Int(f"charge_{atom.GetIdx()}_{possible_charge}")
            charge_variables[(atom.GetIdx(), possible_charge)] = symbol
            # Booleanity guarded/relaxed
            aidx = atom.GetIdx()
            if aidx in excluded_atoms:
                c_bool = BoolVal(True)
                c_bool_name = f"bool_charge_relaxed_{aidx}_{possible_charge}"
            else:
                base_c_bool = Or(symbol == 0, symbol == 1)
                if enable_atom_activation:
                    c_bool = Implies(atom_activation_vars[aidx], base_c_bool)
                    c_bool_name = f"bool_charge_{aidx}_{possible_charge}_conditional"
                else:
                    c_bool = base_c_bool
                    c_bool_name = f"bool_charge_{aidx}_{possible_charge}"
            charge_is_bool.append(c_bool)
            named_constraints[c_bool_name] = c_bool
            penalty_expressions.append(np.abs(possible_charge) * symbol)
        # charge exclusivity
        atom_idx = atom.GetIdx()
        
        if atom_idx in excluded_atoms:
            # For excluded atoms, just require non-negative values (relaxed constraint)
            relaxed_constraint = And([charge_variables[(atom_idx, possible_charge)] >= 0 for possible_charge in possible_charges[atom.GetSymbol()]])
            constraint = relaxed_constraint
            constraint_name = f"charge_relaxed_{atom_idx}"
        else:
            # For non-excluded atoms, enforce strict exclusivity
            base_constraint = Sum([charge_variables[(atom_idx, possible_charge)] for possible_charge in possible_charges[atom.GetSymbol()]]) == 1
            
            if enable_atom_activation:
                # Only enforce this constraint if the atom is active
                constraint = Implies(atom_activation_vars[atom_idx], base_constraint)
                constraint_name = f"charge_exclusivity_{atom_idx}_conditional"
            else:
                constraint = base_constraint
                constraint_name = f"charge_exclusivity_{atom_idx}"
        
        charge_exclusivity.append(constraint)
        named_constraints[constraint_name] = constraint

    # variable for each possible geometry for each atom
    geometry_variables = {}
    geometry_implications = []
    geometry_exclusivity = []
    geometry_is_bool = []
    allowed_bond_constraints = []
    for atom in considered_atoms:
        for (neighbor_count, atom_charge), angle_to_info in possible_geometries[atom.GetSymbol()].items():
            for possible_angles, info in angle_to_info.items():
                geometry_variables[(atom.GetIdx(), neighbor_count, atom_charge, tuple(possible_angles))] = Int(f"geometry_{atom.GetIdx()}_{neighbor_count}_{atom_charge}_{tuple(possible_angles)}")
                # Geometry booleanity guarded/relaxed
                gvar = geometry_variables[(atom.GetIdx(), neighbor_count, atom_charge, tuple(possible_angles))]
                gaidx = atom.GetIdx()
                if gaidx in excluded_atoms:
                    g_bool = BoolVal(True)
                    g_bool_name = f"bool_geom_relaxed_{gaidx}_{neighbor_count}_{atom_charge}_{tuple(possible_angles)}"
                else:
                    base_g_bool = Or(gvar == 0, gvar == 1)
                    if enable_atom_activation:
                        g_bool = Implies(atom_activation_vars[gaidx], base_g_bool)
                        g_bool_name = f"bool_geom_{gaidx}_{neighbor_count}_{atom_charge}_{tuple(possible_angles)}_conditional"
                    else:
                        g_bool = base_g_bool
                        g_bool_name = f"bool_geom_{gaidx}_{neighbor_count}_{atom_charge}_{tuple(possible_angles)}"
                geometry_is_bool.append(g_bool)
                named_constraints[g_bool_name] = g_bool

                # determine forbidden angles
                forbidden_angle_constraints = []
                for neighbor1, neighbor2 in combinations([neighbor.GetIdx() for neighbor in atom.GetNeighbors()], 2):
                    if (neighbor1, neighbor2) in all_pairwise_bond_angles[atom.GetIdx()]:
                        angle_deviation = np.absolute(np.array(possible_angles) - all_pairwise_bond_angles[atom.GetIdx()][(neighbor1, neighbor2)])
                        if len(angle_deviation) == 0: angle_deviation = np.array([180])
                        used_deviation = np.min(angle_deviation)
                        if used_deviation > angle_tolerance:
                            forbidden_angle_constraints.append(Or(
                                bond_order_variables[(atom.GetIdx(), neighbor1, 0)] == 1,
                                bond_order_variables[(atom.GetIdx(), neighbor2, 0)] == 1
                            ))
                        else:
                            # we weight based on the angle deviation
                            weighted_deviation = used_deviation * angle_deviation_weighting
                            RealVal(weighted_deviation)
                            total_penalty_expression = Sum([(RealVal(weighted_deviation) * bond_order_variables[(atom.GetIdx(), neighbor1, i)]) for i in range(1, 4)])
                            penalty_expressions.append(Sum([(RealVal(weighted_deviation) * bond_order_variables[(atom.GetIdx(), neighbor1, i)]) for i in range(0, 4)]))

                # bond order/charge/angle implications, guarded/relaxed
                center_idx = atom.GetIdx()
                original_impl = Implies(
                    geometry_variables[(center_idx, neighbor_count, atom_charge, tuple(possible_angles))] == 1,
                    And([
                        Or([
                            And([
                                Sum([bond_order_variables[(center_idx, neighbor.GetIdx(), bond_order)] for neighbor in atom.GetNeighbors()]) == bond_count
                            for bond_order, bond_count in bond_order_info.items()])
                        for bond_order_info in info["bond_orders"]]),
                        charge_variables[(center_idx, atom_charge)] == 1,
                        *forbidden_angle_constraints
                    ])
                )
                if center_idx in excluded_atoms:
                    impl_constraint = BoolVal(True)
                    impl_name = f"geom_impl_relaxed_{center_idx}_{neighbor_count}_{atom_charge}_{tuple(possible_angles)}"
                else:
                    if enable_atom_activation:
                        impl_constraint = Implies(atom_activation_vars[center_idx], original_impl)
                        impl_name = f"geom_impl_{center_idx}_{neighbor_count}_{atom_charge}_{tuple(possible_angles)}_conditional"
                    else:
                        impl_constraint = original_impl
                        impl_name = f"geom_impl_{center_idx}_{neighbor_count}_{atom_charge}_{tuple(possible_angles)}"
                geometry_implications.append(impl_constraint)
                named_constraints[impl_name] = impl_constraint
        # geometry exclusivity
        atom_idx = atom.GetIdx()
        base_geometry_excl = Sum([geometry_variables[(atom_idx, neighbor_count, atom_charge, tuple(possible_angles))]
                                   for (neighbor_count, atom_charge), angle_to_info in possible_geometries[atom.GetSymbol()].items()
                                   for possible_angles, info in angle_to_info.items()]) == 1
        # Relax geometry exclusivity for excluded atoms; otherwise condition on activation if enabled
        if atom_idx in excluded_atoms:
            constraint = BoolVal(True)
            constraint_name = f"geometry_relaxed_{atom_idx}"
        else:
            if enable_atom_activation:
                constraint = Implies(atom_activation_vars[atom_idx], base_geometry_excl)
                constraint_name = f"geometry_exclusivity_{atom_idx}_conditional"
            else:
                constraint = base_geometry_excl
                constraint_name = f"geometry_exclusivity_{atom_idx}"
        geometry_exclusivity.append(constraint)
        named_constraints[constraint_name] = constraint

    # forbid strong angle violations

    formula = And(bond_order_exclusivity + charge_exclusivity + geometry_implications + geometry_exclusivity + bond_order_constraints + geometry_is_bool + charge_is_bool + bond_order_is_bool)
    objective = Sum(penalty_expressions)
    variable_dicts = {
        **bond_order_variables,
        **charge_variables,
        **geometry_variables
    }
    if enable_atom_activation:
        return formula, objective, variable_dicts, named_constraints, atom_activation_vars
    else:
        return formula, objective, variable_dicts, named_constraints

def _resolve_unsat_cores_iteratively(mol, indices, all_pairwise_bond_angles, all_bond_lengths):
    """
    Iteratively resolve unsat cores by softening constraints for problematic atoms until SAT is achieved.
    Returns the set of atoms that need softened constraints to make the system SAT.
    """
    from z3 import Solver, Optimize, sat, unsat
    
    all_atoms = set(indices)
    softened_atoms = set()
    
    import os
    
    while True:
        if os.environ.get('XB_DEBUG_SAT') == '1':
            print(f"\n--- Iteration with {len(softened_atoms)} softened atoms ---")
        
        # Build formula with current softened atoms
        formula, objective, variable_dict, named_constraints, atom_activation_vars = _construct_SMT_formulation(
            mol, list(all_atoms), all_pairwise_bond_angles, all_bond_lengths, 
            enable_atom_activation=True, excluded_atoms=softened_atoms
        )
        
        # Assumption-based SAT check so unsat core contains only atom activations
        s = Solver()
        s.add(formula)
        assumptions = []
        for atom_idx in all_atoms:
            if atom_idx not in softened_atoms and atom_idx in atom_activation_vars:
                assumptions.append(atom_activation_vars[atom_idx])
        
        result = s.check(assumptions)
        if os.environ.get('XB_DEBUG_SAT') == '1':
            print(f"SAT check result: {result}")
        
        if result == sat:
            if os.environ.get('XB_DEBUG_SAT') == '1':
                print("✓ System is now SAT! All constraints are compatible.")
            break
        else:
            # Extract problematic atoms from unsat core (subset of assumptions)
            unsat_core = s.unsat_core()
            if os.environ.get('XB_DEBUG_SAT') == '1':
                print(f"Unsat core: {unsat_core}")
            
            problematic_atoms = set()
            for lit in unsat_core:
                lit_name = str(lit)
                if lit_name.startswith("atom_active_"):
                    atom_idx = int(lit_name.split("_")[2])
                    problematic_atoms.add(atom_idx)
            
            if problematic_atoms:
                atom_symbols = [mol.GetAtomWithIdx(idx).GetSymbol() + str(idx) for idx in sorted(problematic_atoms)]
                if os.environ.get('XB_DEBUG_SAT') == '1':
                    print(f"Softening constraints for problematic atoms: {atom_symbols}")
                
                # Add problematic atoms to softened set
                softened_atoms.update(problematic_atoms)
            else:
                if os.environ.get('XB_DEBUG_SAT') == '1':
                    print("No atom activation constraints in unsat core - this means the constraints themselves are contradictory")
                    print("This usually happens when the system is fundamentally unsolvable even with relaxed constraints")
                    print("Breaking out of constraint softening loop")
                break
    
    if os.environ.get('XB_DEBUG_SAT') == '1':
        print(f"\n=== CONSTRAINT SOFTENING COMPLETE ===")
        print(f"Softened atoms: {len(softened_atoms)}")
        if softened_atoms:
            softened_symbols = [mol.GetAtomWithIdx(idx).GetSymbol() + str(idx) for idx in sorted(softened_atoms)]
            print(f"Softened atom symbols: {softened_symbols}")
        print(f"Strict atoms: {len(all_atoms - softened_atoms)}")
    
    return softened_atoms, all_atoms - softened_atoms

def _get_connected_subgraphs(mol, indices):
    """
    Given an RDKit molecule and a set of atom indices, identifies the connected subgraphs
    in the subgraph induced by the indices.

    Args:
        mol (rdkit.Chem.rdchem.Mol): The molecule.
        indices (list or set of int): The atom indices to consider.

    Returns:
        list of sets of int: Each set is a connected component of atom indices.
    """
    # Ensure indices are in a set for faster lookup
    indices_set = set(indices)
    
    # Build adjacency list for the induced subgraph
    adjacency = {idx: [] for idx in indices_set}
    for bond in mol.GetBonds():
        idx1 = bond.GetBeginAtomIdx()
        idx2 = bond.GetEndAtomIdx()
        # Include the bond if both atoms are in the specified indices
        if idx1 in indices_set and idx2 in indices_set:
            adjacency[idx1].append(idx2)
            adjacency[idx2].append(idx1)

    # Initialize variables for DFS
    visited = set()
    connected_components = []

    # Perform DFS to find connected components
    for idx in indices_set:
        if idx not in visited:
            stack = [idx]
            component = set()
            while stack:
                current_idx = stack.pop()
                if current_idx not in visited:
                    visited.add(current_idx)
                    component.add(current_idx)
                    # Add neighbors to stack
                    for neighbor in adjacency[current_idx]:
                        if neighbor not in visited:
                            stack.append(neighbor)
            connected_components.append(component)
    return connected_components

def _get_all_pairwise_bond_angles(mol):
    positions = position(mol)
    pairwise_angles = {}
    for atom in mol.GetAtoms():
        center_atom_angles = {}
        for neighbor1, neighbor2 in combinations([neighbor.GetIdx() for neighbor in atom.GetNeighbors()], 2):
            angle = calculate_angle(positions[neighbor1], positions[atom.GetIdx()], positions[neighbor2])
            center_atom_angles[(neighbor1, neighbor2)] = angle
            center_atom_angles[(neighbor2, neighbor1)] = angle
        pairwise_angles[atom.GetIdx()] = center_atom_angles
    return pairwise_angles

def _get_violating_atom_clusters(mol, all_pairwise_bond_angles, angle_tolerance = 10, add_hydrogens = False):
    positions = position(mol)

    violating_atoms = set()
    considerable_geometries = {}
    # id violating atoms
    for atom in mol.GetAtoms():
        pairwise_bond_angles = all_pairwise_bond_angles[atom.GetIdx()]
        found_angles = set(pairwise_bond_angles.values())

        # get geometry
        this_possible_geometries = possible_geometries.get(atom.GetSymbol(), {}).get((len(atom.GetNeighbors()), int(atom.GetFormalCharge())), {})
        if len(possible_geometries) == 0:
            violating_atoms.add(atom.GetIdx())
            continue

        geometry_violation = True
        #print(f"Checking atom {atom.GetSymbol()}{atom.GetIdx()} with {len(atom.GetNeighbors())} neighbors and charge {atom.GetFormalCharge()} and possible_geometries {this_possible_geometries}")
        for possible_angles, possible_geometry_info in this_possible_geometries.items():
            # check if all of the angles are within the allowed range
            if not ((len(found_angles) == 0) and (len(possible_angles) == 0)):
                try:
                    angle_deviations = np.max(np.absolute(np.fromiter(found_angles, float, len(found_angles))[None, :] - np.array(possible_angles)[:, None]))
                except ValueError:
                    # if the found angles are not in the possible angles, this geometry is not possible
                    #print(f"Failed for atom {atom.GetSymbol()}{atom.GetIdx()} and possible angles {possible_angles}: as not all angles {found_angles} are in {possible_angles}")
                    continue
                if angle_deviations > angle_tolerance:
                    # if not, this geometry is not possible
                    #print(f"Failed for atom {atom.GetSymbol()}{atom.GetIdx()} and possible angles {possible_angles}: as not all angles {found_angles} are within the allowed range")
                    continue

            # check bond orders
            bond_order_violation = True
            for bond_order in possible_geometry_info.get("bond_orders", []):
                bond_order = dict(bond_order)
                for bond in atom.GetBonds():
                    try:
                        #print(bond.GetBondTypeAsDouble())
                        bond_order[bond.GetBondTypeAsDouble()] -= 1
                    except KeyError:
                        #print(f"Failed for atom {atom.GetSymbol()}{atom.GetIdx()} with bond {bond.GetIdx()} as {bond.GetBondTypeAsDouble()} is not in {bond_order}")
                        continue
                if not any(bond_order.values()):
                    bond_order_violation = False
                    break
                #print(f"Failed for atom {atom.GetSymbol()}{atom.GetIdx()} with bond {bond.GetIdx()} as bond order violation is {bond_order}")
            if bond_order_violation:
                #print(f"Failed for atom {atom.GetSymbol()}{atom.GetIdx()} and possible angles {possible_angles}: as no valid bond order was found")
                continue
            # if we reach this point, the geometry is possible
            geometry_violation = False
            considerable_geometries[atom.GetIdx()] = possible_angles
        
        if geometry_violation:
            #print(f"Failed for atom {atom.GetSymbol()}{atom.GetIdx()}")
            violating_atoms.add(atom.GetIdx())
    return violating_atoms