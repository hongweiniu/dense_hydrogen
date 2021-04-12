'''
笔记
----------
计算相关量
'''
import numpy as np


def cal_rs(atoms):
    '''
    功能
    ----------
    计算atoms(单帧)的rs

    参数
    ----------
    atoms: ase的atoms对象(单帧)

    返回值
    ----------
    rs: float
    '''
    system_size = len(atoms.get_positions())
    volume = atoms.get_volume()
    density = system_size/volume
    rs = (3./(4.*np.pi*density))**(1./3.)/0.529
    return rs
