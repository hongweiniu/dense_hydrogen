'''
笔记
----------
和layer有关的一些函数
'''
import numpy as np


def plot_layer_stacking(atoms_com,num_of_layers,output_dir,height=5):
    '''
    功能
    ----------
    画com的堆叠情况

    参数
    ----------
    atoms_com: ASE中的atoms对象(单帧)
    num_of_layers: 层数
    output_dir: 输出文件夹
    height: z坐标直方图峰的最小值

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
    myfun.mkdir("%s/layer_pos"%output_dir)
    num_bins = 300
    delta = 1
    layer_pos = list()
    atoms_com.wrap()
    positions = atoms_com.get_positions()
    cell = atoms_com.get_cell()
    cell_z = cell[2][2]
    fig, axe = plt.subplots(1, 1)
    axe.set_xlabel('z (A)')
    axe.set_ylabel('Counts')
    axe.hist(positions[:, 2], bins=num_bins)
    fig.tight_layout()
    fig.savefig("%s/z_hist.png"%output_dir,dpi=320)
    plt.close("all")
    hist, bins = np.histogram(positions[:, 2], bins=np.arange(-delta, cell_z+delta, cell_z/num_bins))
    # print(hist)
    peaks = signal.find_peaks(hist, height=height, distance=num_bins/num_of_layers/1.5)[0]
    # print(peaks)
    if len(peaks)!=num_of_layers:
        print(len(peaks),output_dir)
    z_centers = [0.5*(bins[j]+bins[j+1]) for j in peaks]
    for j in range(len(peaks)):
        layer_confs = extract_layer_confs(atoms_com, z_centers[j]-cell_z/num_of_layers/2, z_centers[j]+cell_z/num_of_layers/2)
        layer_pos.append(layer_confs)
    for j in range(len(peaks)-2):
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
        fig.savefig("%s/layer_pos/layer_pos_%d.png"%(output_dir,j),dpi=320)
        plt.close("all")
