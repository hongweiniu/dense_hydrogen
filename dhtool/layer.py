'''
笔记
----------
和layer有关的一些函数
'''
import numpy as np


def plot_layer_stacking(atoms_com,num_of_layers):
    '''
    功能
    ----------
    画com的堆叠情况

    参数
    ----------
    atoms_com: ASE中的atoms对象(单帧)
    num_of_layers: 层数

    返回值
    ----------
    无
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    from mytool import myfun
    from scipy import signal
    def extract_layer_confs(atoms, zmin, zmax):
        pos = atoms.get_positions()
        z = pos[:, 2]
        sel = (zmin < z) & (z < zmax)
        return pos[sel]
    myfun.mkdir("layer_pos")
    num_bins = 300
    height = 10
    layer_pos = [np.zeros((0, 3))]*num_of_layers
    for i in range(len(atoms_com)):
        positions = atoms_com[i].get_positions()
        cell = atoms_com[i].get_cell()
        cell_z = cell[2][2]
        positions[:, 2] += -np.min(positions[:, 2])+cell_z/num_of_layers/2
        fig, axe = plt.subplots(1, 1)
        axe.set_xlabel('z (A)')
        axe.set_ylabel('Counts')
        axe.hist(positions[:, 2], bins=300)
        fig.tight_layout()
        fig.savefig("zhist.png",dpi=320)
        plt.close("all")
        atoms_com[i].set_positions(positions)
        hist, bins = np.histogram(positions[:, 2], bins=np.arange(0, cell_z+0.001, cell_z/num_bins))
        peaks = signal.find_peaks(hist, height=height, distance=num_bins/num_of_layers/1.5)[0]
        print(len(peaks))
        z_centers = [0.5*(bins[j]+bins[j+1]) for j in peaks]
        if len(peaks) == num_of_layers:
            for j in range(num_of_layers):
                layer_confs = extract_layer_confs(atoms_com[i], z_centers[j]-cell_z/num_of_layers/2, z_centers[j]+cell_z/num_of_layers/2)
                layer_pos[j] = np.r_[layer_pos[j], layer_confs]
            for j in range(num_of_layers-2):
                fig, axe = plt.subplots(1, 1)
                axe.plot(layer_pos[j][:, 0], layer_pos[j][:, 1], linestyle="", marker="_", label="layer_%d" % (j+1))
                axe.plot(layer_pos[j+1][:, 0], layer_pos[j+1][:, 1], linestyle="", marker="|", label="layer_%d" % (j+2))
                axe.plot(layer_pos[j+2][:, 0], layer_pos[j+2][:, 1], linestyle="", marker="x", label="layer_%d" % (j+3))
                axe.set_xlabel('x (A)')
                axe.set_ylabel('y (A)')
                axe.set_xlim(min(positions[:,0])-2,max(positions[:,1])+2)
                axe.set_ylim(min(positions[:,0])-2,max(positions[:,1])+2)
                axe.set_aspect(1.0)
                axe.legend()
                fig.tight_layout()
                fig.savefig("layer_pos/layer_pos_%d.png" % j, dpi=320)
                plt.close("all")
