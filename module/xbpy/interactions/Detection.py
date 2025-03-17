
from . import AnionReceptor, AnionLigand, \
        CationReceptor, CationLigand, \
        HydrogenBondAcceptorReceptor, HydrogenBondDonorLigand, \
        HydrogenBondDonorReceptor, HydrogenBondAcceptorLigand, \
        AromaticProximityLigand, AromaticProximityReceptor
from pathlib import Path

hbond_color = "#4270D4"#"#205F86"
bad_hbond_color = "#4F87FF"#"#D0ECFF"
pi_pi_color = "#000000"#"#6B894E"
cation_dipole_color = "#3F48CC"#"#FF4C0E"
anion_dipole_color = "#3F48CC"#"#0A63FF"
salt_color = "#FF0000"#"#F7D8BD"

class Interactions():

    def __init__(self, receptor):
        self.receptor = receptor
        self.cation_receptor = CationReceptor(receptor)
        self.anion_receptor = AnionReceptor(receptor)
        self.hydrogen_bond_acceptor_receptor = HydrogenBondAcceptorReceptor(receptor)
        self.bad_hydrogen_bond_acceptor_receptor = HydrogenBondAcceptorReceptor(receptor, angle_threshold=60, distance_threshold=3.5)
        self.hydrogen_bond_donor_receptor = HydrogenBondDonorReceptor(receptor)
        self.bad_hydrogen_bond_donor_receptor = HydrogenBondDonorReceptor(receptor, angle_threshold=60, distance_threshold=3.5)
        self.pi_stacking_receptor = AromaticProximityReceptor(receptor)

    def detect_interactions(self, ligand):
        cation_ligand = CationLigand(ligand)
        anion_ligand = AnionLigand(ligand)
        hydrogen_bond_acceptor_ligand = HydrogenBondAcceptorLigand(ligand)
        hydrogen_bond_donor_ligand = HydrogenBondDonorLigand(ligand)
        pi_stacking_ligand = AromaticProximityLigand(ligand)

        # hydrogen_bond_donor_interactions
        hydrogen_bond_donor_interactions = self.hydrogen_bond_donor_receptor.detect_interactions(hydrogen_bond_acceptor_ligand)
        bad_hydrogen_bond_donor_interactions = self.bad_hydrogen_bond_donor_receptor.detect_interactions(hydrogen_bond_acceptor_ligand)
        # remove good interactions from bad interactions
        good_indices = set([(*interaction[:2],) for interaction in hydrogen_bond_donor_interactions])
        bad_hydrogen_bond_donor_interactions = [interaction for interaction in bad_hydrogen_bond_donor_interactions if (*interaction[:2],) not in good_indices]
        # hydrogen_bond_acceptor_interactions
        hydrogen_bond_acceptor_interactions = self.hydrogen_bond_acceptor_receptor.detect_interactions(hydrogen_bond_donor_ligand)
        bad_hydrogen_bond_acceptor_interactions = self.bad_hydrogen_bond_acceptor_receptor.detect_interactions(hydrogen_bond_donor_ligand)
        # remove good interactions from bad interactions
        good_indices = set([(*interaction[:2],) for interaction in hydrogen_bond_acceptor_interactions])
        bad_hydrogen_bond_acceptor_interactions = [interaction for interaction in bad_hydrogen_bond_acceptor_interactions if (*interaction[:2],) not in good_indices]
        # cation-dipole interactions
        cation_dipole_interactions = self.cation_receptor.detect_interactions(hydrogen_bond_acceptor_ligand, distance_threshold=4)
        dipole_cation_interactions = self.hydrogen_bond_acceptor_receptor.detect_interactions(cation_ligand, distance_threshold=4)

        # anion-dipole interactions
        anion_dipole_interactions = self.anion_receptor.detect_interactions(hydrogen_bond_donor_ligand, distance_threshold=4)
        dipole_anion_interactions = self.hydrogen_bond_donor_receptor.detect_interactions(anion_ligand, distance_threshold=4)

        # anion-cation interactions
        anion_cation_interactions = self.anion_receptor.detect_interactions(cation_ligand)
        cation_anion_interactions = self.cation_receptor.detect_interactions(anion_ligand)

        # pi stacking interactions
        pi_stacking_interactions = self.pi_stacking_receptor.detect_interactions(pi_stacking_ligand)

        return {
            "hydrogen_bond_donor_interactions": hydrogen_bond_donor_interactions,
            "hydrogen_bond_acceptor_interactions": hydrogen_bond_acceptor_interactions,
            "bad_hydrogen_bond_donor_interactions": bad_hydrogen_bond_donor_interactions,
            "bad_hydrogen_bond_acceptor_interactions": bad_hydrogen_bond_acceptor_interactions,
            "cation_dipole_interactions": cation_dipole_interactions,
            "dipole_cation_interactions": dipole_cation_interactions,
            "anion_dipole_interactions": anion_dipole_interactions,
            "dipole_anion_interactions": dipole_anion_interactions,
            "anion_cation_interactions": anion_cation_interactions,
            "cation_anion_interactions": cation_anion_interactions,
            "pi_stacking_interactions": pi_stacking_interactions
        }

    def write_interaction_display(self, ligand, filename, prefix = None, with_self = False):
        self.create_interaction_display(ligand, prefix, with_self).write(filename)

    def create_interaction_display(self, ligand, prefix = None, with_self = False, interaction_points = False, return_interactions = False):
        import pymolviz as pmv
        if prefix is None:
            prefix = Path(filename).stem
        interactions = self.detect_interactions(ligand)

        additional_interactions = []
        if with_self:
            ligand_interaction = Interactions(ligand)
            additional_interactions = ligand_interaction.detect_interactions(ligand)
        all_visuals = []

        hbond_points = [[interaction[2], interaction[3]] for interaction in interactions["hydrogen_bond_donor_interactions"]]
        hbond_points +=[[interaction[2], interaction[3]] for interaction in interactions["hydrogen_bond_acceptor_interactions"]]
        if with_self:
            hbond_points +=[[interaction[2], interaction[3]] for interaction in additional_interactions["hydrogen_bond_donor_interactions"]]
            hbond_points +=[[interaction[2], interaction[3]] for interaction in additional_interactions["hydrogen_bond_acceptor_interactions"]]
        if len(hbond_points) > 0:
            all_visuals.append(pmv.Lines(hbond_points, color=hbond_color, name=prefix + " Hydrogen Bonds").as_dotted())
            if interaction_points:
                # generate seperate visuals for receptor and ligand
                all_visuals.append(pmv.Points([interaction[2] for interaction in interactions["hydrogen_bond_donor_interactions"]], color=hbond_color, transparency= 0.9, name=prefix + " Hydrogen Bond Donors Receptor"))
                all_visuals.append(pmv.Points([interaction[2] for interaction in interactions["hydrogen_bond_acceptor_interactions"]], color=hbond_color, transparency= 0.9, name=prefix + " Hydrogen Bond Acceptors Receptor"))
                all_visuals.append(pmv.Points([interaction[3] for interaction in interactions["hydrogen_bond_donor_interactions"]], color=hbond_color, transparency= 0.9, name=prefix + " Hydrogen Bond Donors Ligand"))
                all_visuals.append(pmv.Points([interaction[3] for interaction in interactions["hydrogen_bond_acceptor_interactions"]], color=hbond_color, transparency= 0.9, name=prefix + " Hydrogen Bond Acceptors Ligand"))



        bad_hbond_points = [[interaction[2], interaction[3]] for interaction in interactions["bad_hydrogen_bond_donor_interactions"]]
        bad_hbond_points +=[[interaction[2], interaction[3]] for interaction in interactions["bad_hydrogen_bond_acceptor_interactions"]]
        if with_self:
            bad_hbond_points +=[[interaction[2], interaction[3]] for interaction in additional_interactions["bad_hydrogen_bond_donor_interactions"]]
            bad_hbond_points +=[[interaction[2], interaction[3]] for interaction in additional_interactions["bad_hydrogen_bond_acceptor_interactions"]]
        if len(bad_hbond_points) > 0:
            all_visuals.append(pmv.Lines(bad_hbond_points, color=bad_hbond_color, name=prefix + " Bad Hydrogen Bonds", linewidth=0.035).as_dotted())
            if interaction_points:
                # generate seperate visuals for receptor and ligand
                all_visuals.append(pmv.Points([interaction[2] for interaction in interactions["bad_hydrogen_bond_donor_interactions"]], color=bad_hbond_color, transparency= 0.9, name=prefix + " Bad Hydrogen Bond Donors Receptor"))
                all_visuals.append(pmv.Points([interaction[2] for interaction in interactions["bad_hydrogen_bond_acceptor_interactions"]], color=bad_hbond_color, transparency= 0.9, name=prefix + " Bad Hydrogen Bond Acceptors Receptor"))
                all_visuals.append(pmv.Points([interaction[3] for interaction in interactions["bad_hydrogen_bond_donor_interactions"]], color=bad_hbond_color, transparency= 0.9, name=prefix + " Bad Hydrogen Bond Donors Ligand"))
                all_visuals.append(pmv.Points([interaction[3] for interaction in interactions["bad_hydrogen_bond_acceptor_interactions"]], color=bad_hbond_color, transparency= 0.9, name=prefix + " Bad Hydrogen Bond Acceptors Ligand"))

        cation_dipole_points = [[interaction[2], interaction[3]] for interaction in interactions["cation_dipole_interactions"]]
        cation_dipole_points +=[[interaction[2], interaction[3]] for interaction in interactions["dipole_cation_interactions"]]
        if with_self:
            cation_dipole_points +=[[interaction[2], interaction[3]] for interaction in additional_interactions["cation_dipole_interactions"]]
            cation_dipole_points +=[[interaction[2], interaction[3]] for interaction in additional_interactions["dipole_cation_interactions"]]
        if len(cation_dipole_points) > 0:
            all_visuals.append(pmv.Lines(cation_dipole_points, color=cation_dipole_color, name=prefix + " Cation Dipole Interactions", linewidth=0.065).as_dotted())
            if interaction_points:
                # generate seperate visuals for receptor and ligand
                all_visuals.append(pmv.Points([interaction[2] for interaction in interactions["cation_dipole_interactions"]], color=cation_dipole_color, transparency= 0.9, name=prefix + " Cation Dipole Acceptors Receptor"))
                all_visuals.append(pmv.Points([interaction[2] for interaction in interactions["dipole_cation_interactions"]], color=cation_dipole_color, transparency= 0.9, name=prefix + " Cation Dipole Donors Receptor"))
                all_visuals.append(pmv.Points([interaction[3] for interaction in interactions["cation_dipole_interactions"]], color=cation_dipole_color, transparency= 0.9, name=prefix + " Cation Dipole Donors Ligand"))
                all_visuals.append(pmv.Points([interaction[3] for interaction in interactions["dipole_cation_interactions"]], color=cation_dipole_color, transparency= 0.9, name=prefix + " Cation Dipole Acceptors Ligand"))

        anion_dipole_points = [[interaction[2], interaction[3]] for interaction in interactions["anion_dipole_interactions"]]
        anion_dipole_points +=[[interaction[2], interaction[3]] for interaction in interactions["dipole_anion_interactions"]]
        if with_self:
            anion_dipole_points +=[[interaction[2], interaction[3]] for interaction in additional_interactions["anion_dipole_interactions"]]
            anion_dipole_points +=[[interaction[2], interaction[3]] for interaction in additional_interactions["dipole_anion_interactions"]]
        if len(anion_dipole_points) > 0:
            all_visuals.append(pmv.Lines(anion_dipole_points, color=anion_dipole_color, name=prefix + " Anion Dipole Interactions", linewidth=0.065).as_dotted())
            if interaction_points:
                # generate seperate visuals for receptor and ligand
                all_visuals.append(pmv.Points([interaction[2] for interaction in interactions["anion_dipole_interactions"]], color=anion_dipole_color, transparency= 0.9, name=prefix + " Anion Dipole Acceptors Receptor"))
                all_visuals.append(pmv.Points([interaction[2] for interaction in interactions["dipole_anion_interactions"]], color=anion_dipole_color, transparency= 0.9, name=prefix + " Anion Dipole Donors Receptor"))
                all_visuals.append(pmv.Points([interaction[3] for interaction in interactions["anion_dipole_interactions"]], color=anion_dipole_color, transparency= 0.9, name=prefix + " Anion Dipole Donors Ligand"))
                all_visuals.append(pmv.Points([interaction[3] for interaction in interactions["dipole_anion_interactions"]], color=anion_dipole_color, transparency= 0.9, name=prefix + " Anion Dipole Acceptors Ligand"))

        anion_cation_points = [[interaction[2], interaction[3]] for interaction in interactions["anion_cation_interactions"]]
        anion_cation_points +=[[interaction[2], interaction[3]] for interaction in interactions["cation_anion_interactions"]]
        if with_self:
            anion_cation_points +=[[interaction[2], interaction[3]] for interaction in additional_interactions["anion_cation_interactions"]]
            anion_cation_points +=[[interaction[2], interaction[3]] for interaction in additional_interactions["cation_anion_interactions"]]
        if len(anion_cation_points) > 0:
            all_visuals.append(pmv.Lines(anion_cation_points, color=salt_color, name=prefix + " Salt Bridges", linewidth=0.08).as_dotted())
            if interaction_points:
                # generate seperate visuals for receptor and ligand
                all_visuals.append(pmv.Points([interaction[2] for interaction in interactions["anion_cation_interactions"]], color=salt_color, transparency= 0.9, name=prefix + " Salt Bridge Anions Receptor"))
                all_visuals.append(pmv.Points([interaction[2] for interaction in interactions["cation_anion_interactions"]], color=salt_color, transparency= 0.9, name=prefix + " Salt Bridge Cations Receptor"))
                all_visuals.append(pmv.Points([interaction[3] for interaction in interactions["anion_cation_interactions"]], color=salt_color, transparency= 0.9, name=prefix + " Salt Bridge Cations Ligand"))
                all_visuals.append(pmv.Points([interaction[3] for interaction in interactions["cation_anion_interactions"]], color=salt_color, transparency= 0.9, name=prefix + " Salt Bridge Anions Ligand"))

        pi_stacking_points = [[interaction[2], interaction[3]] for interaction in interactions["pi_stacking_interactions"]]
        if with_self:
            pi_stacking_points +=[[interaction[2], interaction[3]] for interaction in additional_interactions["pi_stacking_interactions"]]
        if len(pi_stacking_points) > 0:
            all_visuals.append(pmv.Lines(pi_stacking_points, color=pi_pi_color, name=prefix + " Pi Stacking Interactions", linewidth=0.1).as_dotted())
            if interaction_points:
                # generate seperate visuals for receptor and ligand
                all_visuals.append(pmv.Points([interaction[2] for interaction in interactions["pi_stacking_interactions"]], color=pi_pi_color, transparency= 0.9, name=prefix + " Pi Stacking Donors Receptor"))
                all_visuals.append(pmv.Points([interaction[3] for interaction in interactions["pi_stacking_interactions"]], color=pi_pi_color, transparency= 0.9, name=prefix + " Pi Stacking Acceptors Ligand"))

        visuals = pmv.Group(all_visuals, name = prefix + " Interactions_Group")
        if return_interactions:
            return visuals, interactions
        return visuals
