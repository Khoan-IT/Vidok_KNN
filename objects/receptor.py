from rdkit.Chem.rdmolfiles import SDMolSupplier
from biopandas.pdb import PandasPdb
from openbabel import openbabel
from objects import construct_hydrogen_bonding_info
import pandas as pd
import os

obConversion = openbabel.OBConversion()

class Receptor():
    def __init__(self, pdb_file):
        self.pdb_file = pdb_file
        self.receptor_df = PandasPdb().read_pdb(self.pdb_file)
        self.receptor_atoms = self.receptor_df.df["ATOM"]
        self.receptor_atoms.sort_values(by=['atom_number'], inplace=True, ignore_index=True)

        fex = self.pdb_file.split(".")[-1]
        self.sdf_file = os.path.dirname(self.pdb_file) + "/SDF/" + os.path.basename(self.pdb_file).split(".")[0] + ".sdf"
        if not os.path.exists(self.sdf_file):
            obConversion.SetInAndOutFormats(fex, "sdf")

            obConversion.OpenInAndOutFiles(self.pdb_file, self.sdf_file)
            obConversion.Convert()
            obConversion.CloseOutFile()

        self.receptor_mol = SDMolSupplier(self.sdf_file)[0]
        self.receptor_bonding_info = construct_hydrogen_bonding_info(self.receptor_mol)

    def get_coords(self):
        coords = self.receptor_atoms[["x_coord", "y_coord", "z_coord"]]
        return coords

    def get_atoms_by_idx(self, list_idx):
        return self.receptor_atoms.iloc[list_idx]

    def get_name(self):
        return os.path.basename(self.pdb_file)

    def save_atoms(self, saving_file):
        self.receptor_atoms.to_csv(saving_file)

    def df_idx_to_rdkit_idx(self, df_idxs):
        return self.receptor_atoms["atom_number"].iloc[df_idxs].to_numpy()
    
    def rdkit_idx_to_df_idx(self, rdkit_idx):
        return self.receptor_atoms[self.receptor_atoms["atom_number"] == rdkit_idx].index.to_numpy()[0]

    @classmethod
    def df_idx_to_atom_num(cls, df_idx):
        return df_idx+1