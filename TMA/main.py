# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 11:32:04 2023

@author: LiHaoJie
"""

from ase import Atoms
from ase.io import read,write
import numpy as np
import pandas as pd
# from ase.build import read_xsd
from ase.visualize import view
from copy import deepcopy

# def translation(atoms, center_atom_index, new_position) :
#     vector = new_position - atoms.positions[center_atom_index]
#     atoms.positions = atoms.positions + vector

def translation(row) :
    new_position = np.array(row['position'])
    # ???
    atoms = deepcopy(row['ligand'][0])
    center_atom_index = row['ligand'][1]
    vector = new_position - atoms.positions[center_atom_index]
    atoms.positions = atoms.positions + vector
    return atoms

def run(file_path, dict_ligand, lattice_vector):
    with open(file_path, 'r', encoding='utf-8') as file:

        # f = file.read()
        flag = 0 #新的时间步
        row = 0 #当前所在行数
        natoms = -1 #总原子数
        time = -1
        data = []
        for i,line in enumerate(file):

            line = line.strip()
            
            if i == 3:
                natoms = int(line)
                continue
        
            if 'ITEM: ATOMS' in line:
                time += 1
                flag = 1
                row = 0
                columns = line.split()[2:]
                continue
            
            if flag == 1:
                str_list = line.split()
                values = [float(v) for v in str_list]
                data.append(values)
                row += 1
            
            if row==natoms:
                df = pd.DataFrame(data)
                df.columns = columns
                flag = 0
                row = 0
                atoms = get_atoms(dict_ligand, df, lattice_vector)
                write('./xyz/%s.cif'%str(time),atoms)
                # del atoms
                del data
                data = []
                
                
def get_atoms(dict_ligand,df,lattice_vector):
    position = {}
    # print(lattice_vector[0])
    # print(df['x'])
    # df['x'] = df['x'] * lattice_vector[0]
    # df['y'] = df['y'] * lattice_vector[1]
    # df['z'] = df['z'] * lattice_vector[2]
    df['position'] = ((df['x'].values[:, np.newaxis] * lattice_vector[0]) + (df['y'].values[:, np.newaxis] * lattice_vector[1]) + (df['z'].values[:, np.newaxis] * lattice_vector[2])).tolist()
    df['ligand'] = df['i1'].map(dict_ligand)
    df = df.dropna()
    df['all_atoms'] = df.apply(translation, axis=1)

    atoms = df['all_atoms'].sum()  
    atoms.set_cell(lattice_vector)
    return atoms
            
            
# atoms = read('./xsd/H2O.xsd')

# view(atoms)
# new_position = np.array([0,0,0])
# translation(atoms, 0, new_position)

# view(atoms)





if __name__ == "__main__":
    # H2O = read('./mol/H2O.mol')
    # OH = read('./mol/OH.mol')
    # ZnEtO = read('./mol/ZnEtO.mol')
    # ZnEtOH = read('./mol/ZnEtOH.mol')
    # ZnEt = read('./mol/ZnEt.mol')
    # ZnEtOH2 = read('./mol/ZnEtOH2.mol')
    # ZnEt2OH = read('./mol/ZnEt2OH.mol')
    # ZnEt2O = read('./mol/ZnEt2O.mol')
    O = Atoms('O', positions=[(0, 0, 0)])
    # Zn = Atoms('Zn', positions=[(0, 0, 0)])
    
    # Hf = read('./mol/Hf.mol')
    TaX5O = read('./mol/TaX5O.mol')
    TaX4O = read('./mol/TaX4O.mol')
    TaX3O = read('./mol/TaX3O.mol')
    TaO = read('./mol/TaO.mol')
    TaX = read('./mol/TaX.mol')
    TaX5OH = read('./mol/TaX5OH.mol')
    O3Ta = read('./mol/O3Ta.mol')
    OH = read('./mol/OH.mol')
    OTa = read('./mol/OTa.mol')
    Ta = Atoms('Ta', positions=[(0, 0, 0)])
    # dict_ligand = {1:(O,0),3:(H2O,0),11:(Zn,0),2:(OH,0),7:(ZnEtO,0),
    #                8:(ZnEtOH,0),12:(ZnEt,0),6:(ZnEtOH2,0),5:(ZnEt2OH,0),4:(ZnEt2O,0)}

    dict_ligand = {1:(O,0),2:(TaX5O,0),3:(TaX4O,0),4:(TaO,0),5:(TaX,0),6:(Ta,0),7:(OTa,0),
                    8:(OH,0),9:(TaX5OH,0),10:(TaX3O,0),11:(O3Ta,0)}
    lattice_vector = np.array([[19.7199001312,0,0],[-9.8599552862,25.9003940686,0],[0,0,63.1218986511]])
    run('./dump.ald',dict_ligand,lattice_vector)







