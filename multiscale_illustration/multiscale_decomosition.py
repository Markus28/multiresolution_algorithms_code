import numpy as np
import matplotlib.pyplot as plt
from denoising import project


def multi_scale_sequence(image, lambdas, tau=1/8):
    ls = 1/(2*lambdas)
    q0 = image - project(image, ls[0], tau)
    partial_sums = [q0]
    summands = [np.copy(q0)]

    for l in ls[1:]:
        summands.append(image - partial_sums[-1] - project(image - partial_sums[-1], l, tau))
        partial_sums.append(partial_sums[-1] + summands[-1])

    return summands, partial_sums



if __name__ == "__main__":
    test_image = 255*plt.imread("ethz_main.png")[:, :, 0]
    a = 1/500
    N = 5
    lambdas = a*3**np.linspace(0, N, num=N+1)
    summands, partial_sums = multi_scale_sequence(test_image, lambdas)

    fig, ax = plt.subplots(ncols = N, nrows = 2)

    
    for i in range(N):
        ax[0, i].axis("off")
        ax[1, i].axis("off")
        ax[1, i].imshow(partial_sums[i], cmap="gray")
        ax[0, i].imshow(summands[i], cmap="gray")
        ax[0, i].set_title(f"$q_{i}$", verticalalignment='top', y=-0.15)
        ax[1, i].set_title(f"""$\\tilde{{q}}_{i}$ \n $\lambda_{i} \\approx {round(lambdas[i], 3)}$""", verticalalignment='top', y=-0.15)

    fig.tight_layout()
    #plt.savefig("multiscale_decomposition.png", dpi=1200)
    plt.savefig("multiscale_decomposition.pdf", dpi=1200, bbox_inches='tight')
    plt.show()
