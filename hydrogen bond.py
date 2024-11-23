import mdtraj as md
import os
import pandas as pd

amino_acids = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
} 

protein = ['T7']
for i,p in enumerate(protein):
    folder_path = f"{p}_pdbs/"
    all_hbond = []
    file_names = [f for f in os.listdir(folder_path) if "relaxed" in f and "unrelaxed" not in f]

    md_WT_protein = md.load(p+'_wt.pdb')
    hbonds_WT = md.baker_hubbard(md_WT_protein, periodic=False)
    label_WT = lambda hbond : '%s -- %s' % (md_WT_protein.topology.atom(hbond[0]), md_WT_protein.topology.atom(hbond[2]))

    for file_name in file_names:
        print(file_name)
        protein_name = os.path.join(folder_path,file_name)
        temp = file_name.split('.')[0].split('_')
        o_res = temp[1]
        site = temp[2]
        m_res = temp[3]
    # Load the protein structure
        md_protein = md.load(protein_name)
        hbonds = md.baker_hubbard(md_protein, periodic=False)
        label = lambda hbond : '%s -- %s' % (md_protein.topology.atom(hbond[0]), md_protein.topology.atom(hbond[2]))
        pair_interest_WT = []
        pair_interest = []

        for hbond in hbonds:
            pair = label(hbond)
            res_interest = amino_acids[m_res]+site+'-'
            if res_interest in pair:
                pair_interest.append(pair)
        for hbond in hbonds_WT:
            pair = label_WT(hbond)
            res_interest = amino_acids[o_res]+site+'-'
            if res_interest in pair:
                pair_interest_WT.append(pair)
        hbond_str = str(len(pair_interest_WT))+'//'+str(len(pair_interest))
        hbond_change = len(pair_interest)-len(pair_interest_WT)
        all_hbond.append([o_res+site+m_res,hbond_str,hbond_change])

    all_df = pd.DataFrame(all_hbond,columns=['mutant','num_hbond(WT//mutant)','change_hbond'])
    all_df = all_df.sort_values(by=['mutant'],ascending=[True])
    all_df.to_csv(f'{p}_hbond.csv',index=False)
