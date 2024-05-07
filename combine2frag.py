#!/python
## author: Zhaoqiang Chen
## contact: 18321885908
from tqdm import tqdm
from glob import glob
from rdkit import Chem
from rdkit import Chem
from rdkit.Chem import AllChem
import os

# if os.path.exists('AB'):
#     raise("please delete or rename AB dir")
# else:
#     os.mkdir("AB")

def combine2frag(Amol,Fr,Bmol,Cs):
    combo = Chem.CombineMols(Amol,Bmol)
    Fr_NEI_ID=get_neiid_bysymbol(combo,Fr)
    Cs_NEI_ID=get_neiid_bysymbol(combo,Cs)
    edcombo = Chem.EditableMol(combo)
    edcombo.AddBond(Fr_NEI_ID,Cs_NEI_ID,order=Chem.rdchem.BondType.SINGLE)

    Fr_ID=get_id_bysymbol(combo,Fr)
    edcombo.RemoveAtom(Fr_ID)
    back = edcombo.GetMol()


    Cs_ID=get_id_bysymbol(back,Cs)

    edcombo=Chem.EditableMol(back)
    edcombo.RemoveAtom(Cs_ID)
    back = edcombo.GetMol()
    smi= Chem.MolToSmiles(back)
    return smi


def get_neiid_bysymbol(combo,symbol):
    for at in combo.GetAtoms():
        if at.GetSymbol()==symbol:
            at_nei=at.GetNeighbors()[0]
            return at_nei.GetIdx()
def get_id_bysymbol(combo,symbol):
    for at in combo.GetAtoms():
        if at.GetSymbol()==symbol:
            return at.GetIdx()