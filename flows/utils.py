from configs import *
import os

def ensure_path(path_name):
    if not os.path.exists(path_name):
        os.mkdir(path_name)

def save_results(result_type=None, **kwargs):
    if not result_type:
        return

    saving_path = os.path.join(results_path, result_type)
    ensure_path(saving_path)

    if result_type == "nearest_atoms":
        for k, v in kwargs.items():
            element_path = os.path.join(saving_path, k)
            ensure_path(element_path)

            if k != "lr_mapping":
                for molecule in v:
                    saving_file = os.path.join(element_path, molecule.get_name() + ".csv")
                    molecule.save_atoms(saving_file)

            else:
                for ligand_idx, per_ligand in enumerate(v):
                    ligand_name = kwargs["ligand_atoms"][ligand_idx].get_name()

                    for receptor_idx, per_receptor in enumerate(per_ligand):
                        receptor_name = kwargs["receptor_atoms"][receptor_idx].get_name()
                        mapping_name = ligand_name + "-" + receptor_name + ".csv"
                        saving_file = os.path.join(element_path, mapping_name)

                        with open(saving_file, "w", encoding="utf-8") as f:
                            for latom, ratoms in per_receptor.items():
                                f.write(str(latom)+",")
                                ratoms = [str(x) for x in ratoms["atom_number"]]
                                f.write(",".join(ratoms) + "\n")

    elif result_type == "hydrogen_bonds":
        element_path = os.path.join(saving_path, "hydrogen_bonds")
        ensure_path(element_path)
        
        for ligand_idx, per_ligand in enumerate(kwargs["hydrogen_bonds"]):
            ligand_name = kwargs["ligand_atoms"][ligand_idx].get_name()
            name = ligand_name.split('.')[0]

            for receptor_idx, per_receptor in enumerate(per_ligand):
                mapping_name = name+".csv"
                saving_file = os.path.join(element_path, mapping_name)

                with open(saving_file, "w", encoding="utf-8") as f:
                    for k in per_receptor:
                        # f.write(str(latom)+",")
                        # ratoms = [str(x) for x in ratoms["atom_number"]]
                        # f.write(",".join(ratoms) + "\n")
                        t = list(k)
                        t[1] = t[1]['atom_number']
                        t = [str(i) for i in t]
                        f.write(",".join(t) + "\n")
                        

def read_candidate_link_results(res_file):
    list_bonds = []
    with open(res_file, "r", encoding="utf-8") as rf:
        lines = rf.read().split("\n")
        for line in lines:
            if line != "":
                list_idx = line.split(",")
                list_idx = [int(x) for x in list_idx if x != ""]
                if len(list_idx) < 2:
                    continue
                for r_idx in list_idx[1:]:
                    list_bonds.append((list_idx[0], r_idx))
                    
    return list_bonds