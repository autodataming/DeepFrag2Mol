#!/usr/bin/env python
import click
from rdkit import Chem
from tqdm import tqdm
from glob import glob
from rdkit import Chem
from rdkit import Chem
from rdkit.Chem import AllChem
from combine2frag import combine2frag



@click.command()
@click.option('--scaffold',default='[*]C1=CC=CC=C1',prompt='input scaffold SMILES:',help='the SMILES of the molecular scaffold')
@click.option('--fragfile',default='deepfrag-scores.csv',prompt='input fragment file path and name',help='the file  of the molecular fragments generated by DeepFrag')



def deepfrag2mol(scaffold,fragfile):

    Fr='Fr'
    Cs='Cs'
    fragA=scaffold
    frag_file=fragfile 
    # print(fragA,frag_file)
    outfilename=frag_file.replace('.csv','_mol.smi')
    basename='deepfragR1_H'
    fragA=fragA.replace('*',Fr)
    # print(fragA)
    Amol=Chem.MolFromSmiles(fragA)
    fh=open(frag_file)
    text=fh.readlines()
    count=0
    total=0
    fw=open(outfilename,'w')
    # print(fragA)
    for line in tqdm(text[1:]):
        name,fragB,score=line.split(',')
        # print(fragA,fragB)
        fragB=fragB.replace('*','[%s]'%Cs)
        # print(fragB)
        Bmol=Chem.MolFromSmiles(fragB)
        total+=1
        # Bmol
        try:
            newsmi=combine2frag(Amol,Fr,Bmol,Cs)
            newname="%s%04d"%(basename,int(name))
            newline="%s %s\n"%(newsmi,newname)
            count+=1
            # print(fragA,fragB,newline,)
            fw.write(newline)
        except:
            print("error:",fragA,fragB)
    fw.close()
    print("Congratulations,You have converted %s total fragments to %s molecules, The result is saved in file %s"%(total,count,outfilename))
        
    
if __name__=='__main__':
    '''
    '''
    helpdoc='''
    DeepFrag2Mol: Converting the molecular fragments generated by DeepFrag into complete molecules.
    Usage: deepFrag2mol_CMDgui.py [OPTIONS]

    Options:
    --scaffold TEXT  the SMILES of the molecular scaffold
    --fragfile TEXT  the file  of the molecular fragments generated by DeepFrag
    --help           Show this message and exit.
    '''
    print(helpdoc)
    deepfrag2mol()
    