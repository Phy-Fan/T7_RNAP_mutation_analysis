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

# Charged residues and their charged atoms
acid_residues = {'D': 'ASP', 'E': 'GLU'}
basic_residues = {'R': 'ARG', 'K': 'LYS', 'H': 'HIS'}
charged_atoms = {
    'ARG': ['NH1', 'NH2'],
    'LYS': ['NZ'],
    'HIS': ['ND1', 'NE2'],
    'ASP': ['OD1', 'OD2'],
    'GLU': ['OE1', 'OE2']
}

def calculate_salt_bridges(site, residue, protein_name, distance_cutoff=4.5):
    try:
        cmd.reinitialize()
        cmd.load(protein_name)
    except Exception:
        return 0

    # Determine residue type and corresponding selection
    if residue in basic_residues:
        target_atoms = '+'.join(charged_atoms[basic_residues[residue]])
        opposite_sel = ' or '.join(
            [f"(resn {res} and name {atoms})" for res, atoms in charged_atoms.items() if res in acid_residues.values()]
        )
    elif residue in acid_residues:
        target_atoms = '+'.join(charged_atoms[acid_residues[residue]])
        opposite_sel = ' or '.join(
            [f"(resn {res} and name {atoms})" for res, atoms in charged_atoms.items() if res in basic_residues.values()]
        )
    else:
        return 0

    try:
        target_sel = f"resi {site} and chain A and name {target_atoms}"
        cmd.select("near_residues", f"(chain A and resi {int(site)-1}) or (chain A and resi {int(site)+1})")
        cmd.select("surrounding", f"byres ({target_sel} around {distance_cutoff})")
        cmd.select("salt_bridges", f"surrounding and ({opposite_sel}) and (not near_residues)")
        num_bridges = cmd.count_atoms("salt_bridges") // 2
        return num_bridges
    except Exception:
        return 0

def process_pdb_files(protein_list):

    for p in protein_list:
        folder_path = f"{p}_pdbs/"
        all_salt_bridges = []

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

            wt_salt_bridges = calculate_salt_bridges(site, o_res, wt_protein_name)
            mutant_salt_bridges = calculate_salt_bridges(site, m_res, protein_name)
            salt_bridge_change = mutant_salt_bridges - wt_salt_bridges

            all_salt_bridges.append([
                res_m,
                wt_salt_bridges,
                mutant_salt_bridges,
                salt_bridge_change
            ])

        columns = ['mutant', 'wt_salt_bridges', 'mutant_salt_bridges', 'change_in_salt_bridges']
        salt_bridges_df = pd.DataFrame(all_salt_bridges, columns=columns)
        salt_bridges_df = salt_bridges_df.sort_values(by=['mutant'], ascending=[True])
        output_file = f'{p}_salt_bridges_analysis.csv'
        salt_bridges_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    protein = ['T7']
    process_pdb_files(protein)
