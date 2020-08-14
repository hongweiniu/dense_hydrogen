'''
笔记
----------
和bond有关的一些函数
'''
import pandas as pd
import numpy as np
from ase import io


def cal_bond_list(f_extxyz, frame, bond_length):
    '''
    功能
    ----------
    计算extxyz文件第frame帧的键列表

    参数
    ----------
    f_extxyz: extxyz文件名
    frame: which frame of extxyz to calculate bond list
    bond_length: 最大键长

    返回值
    ----------
    无
    '''
    atoms = io.read(f_extxyz, index=':', format='extxyz')
    print(atoms)

