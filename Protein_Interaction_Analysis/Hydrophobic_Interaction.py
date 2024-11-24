from pymol import cmd
import os
import pandas as pd

# Amino acid single-letter to three-letter code mapping
amino_acids = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}

# Hydrophobic residues and their hydrophobic atoms
hydrophobic_atoms = {
    'ALA': ['CB'],
    'VAL': ['CG1', 'CG2'],
    'LEU': ['CD1', 'CD2'],
    'ILE': ['CG1', 'CG2', 'CD1'],
    'MET': ['CG', 'SD', 'CE'],
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TRP': ['CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    'PRO': ['CB', 'CG', 'CD']
}

def calculate_hydrophobic_interactions(site, residue, protein_name, distance_cutoff=6.5):
    
    try:
        cmd.reinitialize()
        cmd.load(protein_name)
    except Exception:
        return 0

    # Check if residue is hydrophobic
    res_name = amino_acids[residue]
    if res_name not in hydrophobic_atoms:
        return 0

    try:
        # Construct selection statements for target residue and surrounding residues
        target_atoms = '+'.join(hydrophobic_atoms[res_name])
        target_sel = f"resi {site} and chain A and name {target_atoms}"

        surrounding_sels = []
        for res, atoms in hydrophobic_atoms.items():
            atoms_sel = '+'.join(atoms)
            surrounding_sels.append(f"(resn {res} and name {atoms_sel})")
        surrounding_sel = ' or '.join(surrounding_sels)

        cmd.select("near_residues", f"(chain A and resi {int(site)-1}) or (chain A and resi {int(site)+1})")
        cmd.select("surrounding", f"byres ({target_sel} around {distance_cutoff})")
        cmd.select("hydrophobic_contacts", f"surrounding and ({surrounding_sel}) and (not near_residues)")

        num_contacts = cmd.count_atoms("hydrophobic_contacts") // 2
        return num_contacts

    except Exception:
        return 0

def process_pdb_files(protein_list):
    """
    Process a list of PDB files and calculate hydrophobic interaction changes.
    """
    for p in protein_list:
        folder_path = f"{p}_pdbs/"
        all_interactions = []

        if not os.path.exists(folder_path):
            continue
            
        file_names = [f for f in os.listdir(folder_path) if "relaxed" in f and "unrelaxed" not in f]
        
        wt_protein_name = f"{p}_wt.pdb"
        if not os.path.exists(wt_protein_name):
            continue

        for file_name in file_names:
            protein_name = os.path.join(folder_path, file_name)
            temp = file_name.split('.')[0].split('_')
            o_res = temp[1]
            site = temp[2]
            m_res = temp[3]
            res_m = o_res + site + m_res

            wt_interactions = calculate_hydrophobic_interactions(site, o_res, wt_protein_name)
            mutant_interactions = calculate_hydrophobic_interactions(site, m_res, protein_name)
            interaction_change = mutant_interactions - wt_interactions

            all_interactions.append([
                res_m,
                wt_interactions,
                mutant_interactions,
                interaction_change
            ])

        # Define simplified columns for output
        columns = ['mutant', 'wt_interactions', 'mutant_interactions', 'change_in_interactions']
        interactions_df = pd.DataFrame(all_interactions, columns=columns)
        interactions_df = interactions_df.sort_values(by=['mutant'], ascending=[True])
        output_file = f'{p}_hydrophobic_interactions_analysis.csv'
        interactions_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    protein = ['T7']
    process_pdb_files(protein)
