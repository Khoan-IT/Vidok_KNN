import numpy as np

def find_candidate_amino_atom(receptor_coords, ligand_coords, model):
    list_result = []
    list_distance = []
    for i, receptor_coord in enumerate(receptor_coords):
        list_result.append([])
        list_distance.append([])
        model.fit(receptor_coord)

        for ligand_coord in ligand_coords:
            distances, nearest_ligand_atom_idx = model.radius_neighbors(ligand_coord)
            list_result[i].append(nearest_ligand_atom_idx)
            list_distance[i].append(distances)

    return list_result, list_distance
