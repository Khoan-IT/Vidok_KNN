from flows.candidate_link import find_candidate_amino_atom
from flows.post_process import build_la_ra_mapping_by_idx, filter_by_chemical_properties
from flows.utils import save_results
from objects.data_builder import DataBuilder
from models.configs import model_configs
from models import create_model
from utils import *
import configs
import time
import sys
import numpy as np
import glob

@record_timestamp
def main(FLAGS=None):
    if FLAGS is None:
        FLAGS, _ = configs.parser.parse_known_args(args=sys.argv[1:])

    DataBuilder().init(configs.data_path)
    list_receptors = DataBuilder().get_atom_data(name = '', category = 'receptor')
    receptor_coords = DataBuilder().get_atom_data(name = '', category = 'receptor', only_coords=True)
    
    count_error = 0

    list_file = glob.glob('data/ligands/*.pdb')
    name_files = [fn.split('/')[-1] for fn in list_file]
    for index, n in enumerate(name_files):
        try:
            print("Index: {}, Name: {}".format(index,n))
            # Build data
            list_ligands = DataBuilder().get_atom_data(name = n)
            # print(list_ligands)
            ligand_coords = DataBuilder().get_atom_data(name = n, only_coords=True)

            # Initialize models
            print("CREATING MODELS...")
            model = create_model(FLAGS.model, **model_configs[FLAGS.model])

            # Run program
            # start_time = time.time()
            print("RUNNING ALGORITHM...")
            nearest_amino_atom_idx, list_distances = find_candidate_amino_atom(receptor_coords, ligand_coords, model)
            # Post-process result
            print("POST-PROCESSING RESULTS...")
            # la_ra_mapping = build_la_ra_mapping_by_idx(list_receptors, list_ligands, nearest_amino_atom_idx)
            hydrogen_bonds = filter_by_chemical_properties(list_receptors, list_ligands, nearest_amino_atom_idx, list_distances)
            # print(time.time() - start_time)
            # Save result
            print("SAVING RESULTS...")
            # save_results("nearest_atoms", receptor_atoms=list_receptors,
            #                             ligand_atoms=list_ligands,
            #                             lr_mapping=la_ra_mapping)

            save_results("hydrogen_bonds", receptor_atoms=list_receptors,
                                        ligand_atoms=list_ligands,
                                        hydrogen_bonds=hydrogen_bonds)
            time.sleep(5)
        except:
            print('Errors file: {}'.format(n))
            count_error+=1
    print(count_error)

if __name__ == '__main__':
    _, runtime = main()
    print("Total runtine: %.03f" % runtime)
