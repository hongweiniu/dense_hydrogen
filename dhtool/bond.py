'''
笔记
----------
和bond有关的一些函数
'''
import numpy as np
from ase import Atoms


def cal_bond_list(atoms, bond_length):
    '''
    功能
    ----------
    计算atoms(单帧)的键列表

    参数
    ----------
    atoms: ASE中的atoms对象(单帧)
    bond_length: 最大键长

    返回值
    ----------
    bond_list: numpy array
    '''
    distances = atoms.get_all_distances(mic=True)
    bond_list_init = list()
    for i in range(len(distances)):
        nearest_atom = np.argsort(distances[i])[1]
        if atoms.get_distance(i, nearest_atom, mic=True) > bond_length:
            continue
        bond_list_init += [[i, nearest_atom]]
    # print(bond_list_init)
    # print(len(bond_list_init))
    bond_list = list()
    for i in range(len(bond_list_init)):
        if [bond_list_init[i][0], bond_list_init[i][1]] in bond_list_init and [bond_list_init[i][1], bond_list_init[i][0]] in bond_list_init and bond_list_init[i][0] < bond_list_init[i][1]:
            bond_list += [[bond_list_init[i][0], bond_list_init[i][1]]]
    # print(bond_list)
    return bond_list


def move_bonded_atoms(atoms, bond_list, flag):
    '''
    功能
    ----------
    将atoms(单帧)中边界处成键的原子移动到盒子上方或者下方

    参数
    ----------
    atoms: ASE中的atoms对象(单帧)
    bond_list: bond list
    flag: 移动到上方或者下方的标志。'u'或者'd'

    返回值
    ----------
    atoms_moved: ASE中的atoms对象(单帧)
    '''
    lz = atoms.get_cell()[2][2]
    positions = atoms.get_positions()
    for i in bond_list:
        if np.abs(positions[i[0]][2] - positions[i[1]][2]) > lz/2.0:
            if flag == 'u':
                if positions[i[0]][2] > positions[i[1]][2]:
                    positions[i[1]][2] += lz
                else:
                    positions[i[0]][2] += lz
            if flag == 'd':
                if positions[i[0]][2] > positions[i[1]][2]:
                    positions[i[0]][2] -= lz
                else:
                    positions[i[1]][2] -= lz
    atoms_moved = Atoms(symbols=atoms.get_chemical_symbols(), positions=positions, cell=atoms.get_cell(), pbc=atoms.get_pbc())
    return atoms_moved
