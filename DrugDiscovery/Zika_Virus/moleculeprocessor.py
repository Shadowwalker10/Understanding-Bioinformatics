from PaDEL_pywrapper import PaDEL
from typing import Tuple, List
from PaDEL_pywrapper.descriptor import PubchemFP
from rdkit import Chem
from rdkit.Chem import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize

class MoleculeProcessor:
    def __init__(self):
        self.fp = PubchemFP()
        ##Salt remover
        self.salt_remover = SaltRemover.SaltRemover()
        # self.normalizer = rdMolStandardize.Normalize()
        # self.reionizer = rdMolStandardize.Reionize()
        self.padel = PaDEL(descriptors = [self.fp], ignore_3D=True)

    def preprocess_molecule(self, mol):
        if mol is None:
            return None
        ## Remove Salt From Molecule
        self.salt_remover.StripMol(mol, dontRemoveEverything=True)
        ##Standardize Nitro groups
        mol = rdMolStandardize.Normalize(mol)
        ## Reionize the molecule
        mol = rdMolStandardize.Reionize(mol)
        return mol

    def process_smile(self, smiles:List[str]):
        mols = [self.preprocess_molecule(Chem.MolFromSmiles(smile)) for smile in smiles if Chem.MolFromSmiles(smile) is not None]
        return self.padel.calculate(mols)