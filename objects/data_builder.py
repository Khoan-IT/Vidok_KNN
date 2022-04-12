from objects.receptor import Receptor
from objects.ligand import Ligand
import os

class DataBuilder():
    @classmethod
    def init(cls, data_dir):
        cls.data_dir = data_dir
        cls.receptor_dir = os.path.join(cls.data_dir, "receptors")
        cls.ligand_dir = os.path.join(cls.data_dir, "ligands")

    @classmethod
    def check_file_type(cls, filename):
        fex = filename.split(".")[-1]
        return fex in ["pdb", "pdbqt"]

    @classmethod
    def get_atom_data(cls, name, category='ligand',only_coords=False):
        if category == 'ligand':
            list_ligands = []
            list_ligands.append(
                Ligand(os.path.join(cls.ligand_dir, name))
            )
            if only_coords:
                list_ligands = list(map(lambda x: x.get_coords(), list_ligands))

            return list_ligands
        else:
            list_receptors = []
            list_receptors.append(
                Receptor(os.path.join(cls.receptor_dir, 'ZINC000000150863_2.pdb'))
            )
            if only_coords:
                list_receptors = list(map(lambda x: x.get_coords(), list_receptors))
            
            return list_receptors
