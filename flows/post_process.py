from objects.ligand import Ligand
from objects.receptor import Receptor

def build_la_ra_mapping_by_idx(list_receptors, list_ligands, nearest_amino_atom_idx):
    mapping_la_ra = []
    for j, ligand in enumerate(list_ligands):
        mapping_la_ra.append([])
        for i, receptor in enumerate(list_receptors):
            mapping_la_ra[j].append([])

            for latom_idx, list_idx_each_ligand in enumerate(nearest_amino_atom_idx[i][j]):
                if len(list_idx_each_ligand)!=0:
                    for k in list_idx_each_ligand:
                        receptor_atom = list_receptors[i].get_atoms_by_idx(k[0])
                        mapping_la_ra[j][i].append((latom_idx, receptor_atom, k[1], k[2]))
    return mapping_la_ra

def check_hydrogen_bonds(list_candidate_bonds, flatten_list_distance, ligand_bonding_info, receptor_bonding_info):
    ligand_aod = {}
    receptor_aod = {}
    hydrogen_bonds = []
    
    ligand_aod["Acceptor"] = list(y[0] for y in filter(lambda x: x[1] == "Acceptor", ligand_bonding_info))
    ligand_aod["Donor"] = list(y[0] for y in filter(lambda x: x[1] == "Donor", ligand_bonding_info))
    
    receptor_aod["Acceptor"] = list(y[0] for y in filter(lambda x: x[1] == "Acceptor", receptor_bonding_info))
    receptor_aod["Donor"] = list(y[0] for y in filter(lambda x: x[1] == "Donor", receptor_bonding_info))

    
    for bond, distance in zip(list_candidate_bonds, flatten_list_distance):
        latm_idx, ratm_idx = bond
        if latm_idx in ligand_aod["Acceptor"] and ratm_idx in receptor_aod["Donor"] and distance < 3.5:
            hydrogen_bonds.append((latm_idx, ratm_idx, 'Acceptor', 'Donor'))
            
        if latm_idx in ligand_aod["Donor"] and ratm_idx in receptor_aod["Acceptor"] and distance < 4:
            hydrogen_bonds.append((latm_idx, ratm_idx, 'Donor', 'Acceptor'))            
    
    return hydrogen_bonds

def convert_df_idx_to_rdkit_idx(ligand, receptor, list_bonds_by_ligand, list_distances):
    list_candidate_bonds_by_rdkit_idx = []
    flatten_list_distance = []

    ligand_pdb_idxs = ligand.df_idx_to_rdkit_idx(list(range(len(list_bonds_by_ligand))))
    
    for ligand_df_idx, receptor_df_idxs in enumerate(list_bonds_by_ligand):
        flatten_list_distance += list_distances[ligand_df_idx].tolist()
        receptor_pdb_idxs = receptor.df_idx_to_rdkit_idx(receptor_df_idxs)

        for receptor_pdb_idx in receptor_pdb_idxs:
            list_candidate_bonds_by_rdkit_idx.append((ligand_pdb_idxs[ligand_df_idx], receptor_pdb_idx))

    return list_candidate_bonds_by_rdkit_idx, flatten_list_distance

def convert_rdkit_idx_to_df_idx(flatten_list_bonds, ligand, receptor):
    list_bonds_by_df_idx = [[] for x in ligand.ligand_atoms.index]

    for latom_rdkit, ratom_rdkit, l_type, r_type in flatten_list_bonds:
        latom_df_idx = ligand.rdkit_idx_to_df_idx(latom_rdkit)
        ratom_df_idx = receptor.rdkit_idx_to_df_idx(ratom_rdkit)
        list_bonds_by_df_idx[latom_df_idx].append((ratom_df_idx, l_type, r_type))

    return list_bonds_by_df_idx

def filter_by_chemical_properties(list_receptors, list_ligands, mapping_la_ra, list_distances):
    list_hydrogen_bonds = []
    for ligand_idx, per_ligand in enumerate(list_ligands):
        list_hydrogen_bonds.append([])
        for receptor_idx, per_receptor in enumerate(list_receptors):
            # print(ligand_idx)
            # print(receptor_idx)
            list_candidate_bonds = mapping_la_ra[ligand_idx][receptor_idx]
            list_candidate_bonds, flatten_list_distance = convert_df_idx_to_rdkit_idx(per_ligand, per_receptor, list_candidate_bonds, list_distances[ligand_idx][receptor_idx])

            hydrogen_bonds = check_hydrogen_bonds(list_candidate_bonds, flatten_list_distance,
                                                per_ligand.ligand_bonding_info, 
                                                per_receptor.receptor_bonding_info)
                                                
            hydrogen_bonds = convert_rdkit_idx_to_df_idx(hydrogen_bonds, per_ligand, per_receptor)
            list_hydrogen_bonds[ligand_idx].append(hydrogen_bonds)

    return build_la_ra_mapping_by_idx(list_receptors, list_ligands, list_hydrogen_bonds)