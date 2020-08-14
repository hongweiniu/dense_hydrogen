'''
笔记
----------
和bond有关的一些函数
'''
import numpy as np


def cal_bond_list(atoms, frame, bond_length):
    '''
    功能
    ----------
    计算atoms第frame帧的键列表

    参数
    ----------
    atoms: ASE中的atoms对象
    frame: which frame of extxyz to calculate bond list
    bond_length: 最大键长

    返回值
    ----------
    无
    '''
    distances = atoms[frame].get_all_distances(mic=True)
    bond_list_init = list()
    for i in range(len(distances)):
        nearest_atom = np.argsort(distances[i])[1]
        if atoms[frame].get_distance(i, nearest_atom, mic=True) > bond_length:
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
