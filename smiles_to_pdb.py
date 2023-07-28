"""
Purpose of this script is to take a SMILES string as an input, initialize a 3D structure,
optimize it with a basic force field, and write a PDB containing the coordinates. The 
structures are really only supposed to be thorough enough for docking where the ligand
is flexible anyway.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

output = 'carbon_monoxide'
smiles = 'CO' 
forcefield = 'mmff94'

def optimize(smiles, forcefield='mmff94'):
    # Convert smiles to a rdkit molecule
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # Generate an initial 3D conformation
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

    # Optimize the geometry with the specified force field; default is mmff94
    if forcefield == 'mmff94':
        FF = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
    elif forcefield == 'uff':
        FF = AllChem.UFFGetMoleculeForceField(mol)
    else:
        print("This script doesn't support that force field.")
        return

    FF.Minimize()
    
    writer = Chem.PDBWriter(f'{output}.pdb')
    writer.write(mol)
    writer.close()

if __name__ == "__main__":
    optimize(smiles, forcefield)
