from .interface import PyMOLInterface
from ..interactions import Interactions
import re

def get_first_parentheses_content(s):
    match = re.search(r'\(([^)]*)\)', s)
    return match.group(1) if match else None

def vizualize_ligand_interactions(selection, distance_threshold = 8, only_enabled = True, update_graphics = True):
    """ Creates and loads a pymolviz created visualization of the interactions of the given ligand with its surroundings in the current pymol session.
    """
    from pymol import cmd
    ligand_interface = PyMOLInterface(selection, full_objects=False)
    receptor_selection = f"byres {selection} around {distance_threshold}"
    if only_enabled:
        receptor_selection += " and enabled"
    receptor_interface = PyMOLInterface(receptor_selection, full_objects=False)
    receptor_interactions = {object_name: Interactions(molecule) for object_name, molecule in receptor_interface.molecules.items()}

    from ..rdutil import write_molecules
    # write receptors
    for name, mol in ligand_interface.molecules.items():
        write_molecules(mol, name + "_ligand.sdf")
    # for each ligand and each receptor we can now compute the interactions
    interaction_visuals = {}
    found_interactions = {}
    for ligand_object_name, ligand_molecule in ligand_interface.molecules.items():
        for receptor_object_name, receptor_interaction in receptor_interactions.items():
            ligand_name = get_first_parentheses_content(ligand_object_name)
            receptor_name = get_first_parentheses_content(receptor_object_name)
            interaction_visuals[(ligand_object_name, receptor_object_name)], found_interactions[(ligand_object_name, receptor_object_name)] = receptor_interaction.create_interaction_display(ligand_molecule, prefix = f"int_{ligand_name}_{receptor_name}", return_interactions = True)

    for interaction_visual in interaction_visuals.values():
        interaction_visual.load()
        interaction_visual.write(interaction_visual.name + ".py")

    # additionally we want to only show the interacting sidechains/bbs
    if update_graphics:
        cmd.hide(f"everything", f"({selection}) or (bymol {receptor_selection})")
        # completely show ligand
        cmd.show(f"sticks", selection)

        #show receptor as cartoon
        cmd.show(f"cartoon", f"bymol {receptor_selection}")

        # and interacting side chains as sticks
        for (ligand_object_name, receptor_object_name), interactions in found_interactions.items():
            #get interacting atom indices
            for interaction_type, this_interactions in interactions.items():
                for interaction in this_interactions:
                    atom_indices = interaction[0]
                    try:
                        atom_indices[0]
                    except:
                        atom_indices = [atom_indices]
                    # convert back to pymol selection
                    pymol_query = receptor_interface.get_pymol_query(receptor_object_name, atom_indices)
                    cmd.pseudoatom("residue_labels", f"{pymol_query}", label=f"{cmd.get_model(pymol_query).atom[0].resn} {cmd.get_model(pymol_query).atom[0].resi}")
                    if cmd.select(f"({pymol_query}) and backbone"):
                        cmd.show(f"sticks", f"(({pymol_query}) extend 2) and backbone")
                    else:
                        cmd.show(f"sticks", f"((byres ({pymol_query})) and sidechain) extend 1")
        cmd.color("wheat", f"polymer and e. C")
        cmd.set("float_labels", 1)
        cmd.set("label_bg_color", "white")
        cmd.set("label_font_id", 5)
        cmd.hide("everything", "e. H and (neighbor e. C)")
        cmd.zoom(selection, animate = -1)