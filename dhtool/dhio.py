'''
笔记
----------
和io有关的一些函数
'''
import numpy as np
from ase import Atoms
from dhtool import bond


def gen_com_atoms(atoms, bond_length):
    '''
    功能
    ----------
    使用h2的atoms(单帧)生成com的atoms(单帧)

    参数
    ----------
    atoms: h2的atoms对象(单帧)
    bond_length: 最大键长

    返回值
    ----------
    atoms_com: com的atoms对象(单帧)
    '''
    bond_list = bond.cal_bond_list(atoms, bond_length)
    positions = atoms.get_positions()
    positions_com = np.zeros([len(bond_list), 3])
    for i in range(len(bond_list)):
        distance_vector = atoms.get_distance(bond_list[i][0], bond_list[i][1], mic=True, vector=True)
        positions_com[i] = positions[bond_list[i][0]] + 0.5*distance_vector
    return Atoms('H%d' % (len(bond_list)), positions=positions_com, cell=atoms.cell, pbc=[1, 1, 1])
    