from rdkit.Chem.Lipinski import *
import os

class _ChemicalFeaturesFactory:
    """This is a singleton class for RDKit base features."""
    _instance = None

    @classmethod
    def get_instance(cls):
        try:
            from rdkit import RDConfig
            from rdkit.Chem import ChemicalFeatures
        except ModuleNotFoundError:
            raise ImportError("This class requires RDKit to be installed.")

        if not cls._instance:
            fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
            cls._instance = ChemicalFeatures.BuildFeatureFactory(fdefName)
        return cls._instance

def construct_hydrogen_bonding_info(mol):
    """Construct hydrogen bonding infos about a molecule.
    Parameters
    ---------
    mol: rdkit.Chem.rdchem.Mol
    RDKit mol object
    Returns
    -------
    List[Tuple[int, str]]
    A list of tuple `(atom_index, hydrogen_bonding_type)`.
    The `hydrogen_bonding_type` value is "Acceptor" or "Donor".
    """
    factory = _ChemicalFeaturesFactory.get_instance()
    feats = factory.GetFeaturesForMol(mol)
    hydrogen_bonding = []
    for f in feats:
        hydrogen_bonding.append((f.GetAtomIds()[0], f.GetFamily()))
    return hydrogen_bonding
