import mdtraj as md
import os
import pandas as pd

# Amino acid single-letter to three-letter code mapping
amino_acids = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}

def calculate_hbonds(site, o_res, m_res, wt_protein_path, mutant_protein_path):

    # Load wild-type and mutant protein structures
    md_WT_protein = md.load(wt_protein_path)
    md_mutant_protein = md.load(mutant_protein_path)

    # Compute H-bonds
    hbonds_WT = md.baker_hubbard(md_WT_protein, periodic=False)
    hbonds_mutant = md.baker_hubbard(md_mutant_protein, periodic=False)

    # Labeling function
    label_WT = lambda hbond: '%s -- %s' % (md_WT_protein.topology.atom(hbond[0]), md_WT_protein.topology.atom(hbond[2]))
    label_mutant = lambda hbond: '%s -- %s' % (md_mutant_protein.topology.atom(hbond[0]), md_mutant_protein.topology.atom(hbond[2]))

    # Identify relevant H-bonds
    res_interest_WT = amino_acids[o_res] + site + '-'
    res_interest_mutant = amino_acids[m_res] + site + '-'

    pair_interest_WT = [label_WT(hbond) for hbond in hbonds_WT if res_interest_WT in label_WT(hbond)]
    pair_interest_mutant = [label_mutant(hbond) for hbond in hbonds_mutant if res_interest_mutant in label_mutant(hbond)]

    # Calculate counts and changes
    hbond_count_WT = len(pair_interest_WT)
    hbond_count_mutant = len(pair_interest_mutant)
    hbond_change = hbond_count_mutant - hbond_count_WT

    return hbond_count_WT, hbond_count_mutant, hbond_change

def process_pdb_files(protein_list):

    for p in protein_list:
        folder_path = f"{p}_pdbs/"
        all_hbond = []

        if not os.path.exists(folder_path):
            continue

        file_names = [f for f in os.listdir(folder_path) if "relaxed" in f and "unrelaxed" not in f]
        wt_protein_path = f"{p}_wt.pdb"

        if not os.path.exists(wt_protein_path):
            continue

        for file_name in file_names:
            mutant_protein_path = os.path.join(folder_path, file_name)
            temp = file_name.split('.')[0].split('_')
            o_res, site, m_res = temp[1], temp[2], temp[3]
            res_m = o_res + site + m_res

            wt_hbond, mutant_hbond, change_hbond = calculate_hbonds(site, o_res, m_res, wt_protein_path, mutant_protein_path)
            all_hbond.append([res_m, wt_hbond, mutant_hbond, change_hbond])

        # Save results to a CSV file
        columns = ['mutant', 'wt_hbond', 'mutant_hbond', 'change_hbond']
        hbond_df = pd.DataFrame(all_hbond, columns=columns)
        hbond_df = hbond_df.sort_values(by=['mutant'], ascending=True)
        output_file = f"{p}_hbond_analysis.csv"
        hbond_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    protein = ['T7']
    process_pdb_files(protein)
