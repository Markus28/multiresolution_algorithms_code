import numpy as np
from scipy.ndimage.filters import laplace
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
import os

def read_grid(file_name):
    with open(file_name) as f:
        N, M = f.readline().split(",")
        data = np.zeros((int(N), int(M)))
        for line in f:
            i, j, num = line.split(",")
            data[int(i), int(j)] = float(num)

    return data


def transform_for_imshow(A):
    return A.T[::-1]


def visualize_forward():
    location = "data/forward"
    grid_files = []
    for file_name in os.listdir(location):
        if os.path.splitext(file_name)[1] == ".txt":
            grid_files.append(os.path.join(location, file_name))
                

    grid_files = sorted(grid_files, key = lambda f: int(f.split('_')[1].split('.')[0]))

    print(f"Found {len(grid_files)} files...")
    print("Plotting...")

    grids = [transform_for_imshow(read_grid(file_name)) for file_name in grid_files]
    grids = [grid-grids[0] for grid in grids]

    mini = min([np.min(grid) for grid in grids])
    maxi = max([np.max(grid) for grid in grids])

    
    fig,ax = plt.subplots(1)
    rect = patches.Rectangle(((3*grids[0].shape[0])/8, (3*grids[0].shape[1])/8), grids[0].shape[0]/4, grids[0].shape[1]/4, linewidth=1,edgecolor='r',facecolor='none')
    dot = patches.Circle(((1+np.cos(2.3))*grids[0].shape[0]/2, (1-np.sin(2.3))*grids[0].shape[0]/2), radius = grids[0].shape[0]/100, facecolor='r')
    ims = []
    h_r = 2.0/(len(grids) - 1)

    for j, grid in enumerate(grids):
        circle = patches.Circle(((1+np.cos(2.3))*grids[0].shape[0]/2, (1-np.sin(2.3))*grids[0].shape[0]/2), radius = j*h_r*(grids[0].shape[0])/2, linewidth=1, edgecolor='b',facecolor='none')
        ims.append([ax.imshow(grid, vmin = mini, vmax = maxi, animated = True, cmap="coolwarm"), ax.add_patch(rect), ax.add_patch(circle), ax.add_patch(dot)])
    
    ax.axis("off")
    ani = animation.ArtistAnimation(fig, ims, interval=5000/len(grids), blit=True,
                                repeat_delay=0)


    ani.save('visualization.mp4')

    plt.title(r'$\phi - \phi_{*}$')
    plt.show()

    #We want to generate a static plot:
    begin_frame = 67
    frame_skip = 22
    num_frames = 4
    fig, ax = plt.subplots(ncols = num_frames)

    frames = [grids[i*frame_skip + begin_frame] for i in range(num_frames)]
    mini_stat = min([np.min(grid) for grid in frames])
    maxi_stat = max([np.max(grid) for grid in frames])
    
    for i, frame in enumerate(frames):
        rect = patches.Rectangle(((3*grids[0].shape[0])/8, (3*grids[0].shape[1])/8), grids[0].shape[0]/4, grids[0].shape[1]/4, linewidth=1,edgecolor='r',facecolor='none')
        dot = patches.Circle(((1+np.cos(2.3))*grids[0].shape[0]/2, (1-np.sin(2.3))*grids[0].shape[0]/2), radius = grids[0].shape[0]/100, facecolor='r')
        circle = patches.Circle(((1+np.cos(2.3))*grids[0].shape[0]/2, (1-np.sin(2.3))*grids[0].shape[0]/2),
                                radius = (i*frame_skip + begin_frame)*h_r*(grids[0].shape[0])/2, linewidth=1, edgecolor='b',facecolor='none')
        
        ax[i].add_patch(circle)
        ax[i].add_patch(rect)
        ax[i].add_patch(dot)
        ax[i].axis("off")
        ax[i].imshow(frame, vmin = mini_stat, vmax = maxi_stat, cmap="coolwarm")

    fig.tight_layout()
    plt.savefig("forward_problems.pdf", dpi=1200, bbox_inches='tight')
    plt.show()

def visualize_psi():
    location = "data/psi.txt"
    psi = transform_for_imshow(read_grid(location))
    laplace_psi = np.zeros_like(psi)
    laplace_psi[100:300, 100:300] = laplace(psi)[100:300, 100:300]
    fig, ax = plt.subplots(ncols = 2)
    ax[1].imshow(laplace_psi, cmap="coolwarm")
    ax[0].imshow(psi, cmap="coolwarm")
    ax[0].axis("off")
    ax[1].axis("off")
    ax[0].set_title("$\\psi$", verticalalignment='top', y=-0.15)
    ax[1].set_title(f"$\\Delta \\psi$""", verticalalignment='top', y=-0.15)
    fig.tight_layout()
    plt.savefig("psi.pdf", dpi = 1200, bbox_inches='tight')
    plt.show()

def visualize_reconstruction():
    location = "data/absorption.txt"
    reconstruction = transform_for_imshow(read_grid(location))
    plt.imshow(reconstruction, cmap="coolwarm")
    plt.axis("off")
    plt.savefig("reconstruction.pdf", dpi = 1200, bbox_inches='tight')
    plt.show()

if __name__=="__main__":
    visualize_reconstruction()
    visualize_psi()
    visualize_forward()
