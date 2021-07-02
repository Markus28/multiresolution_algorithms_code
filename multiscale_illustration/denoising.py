import numpy as np
import matplotlib.pyplot as plt

def grad(A):
    result = np.zeros((*A.shape, 2))
    result[:,:,0] = np.roll(A, -1, axis=0) - A
    result[:,:,1] = np.roll(A, -1, axis=1) - A
    result[-1, :, 0] = 0
    result[:, -1, 1] = 0
    return result

def div(A):
    X = A[:, :, 0] - np.roll(A[:, :, 0], 1, axis=0)
    Y = A[:, :, 1] - np.roll(A[:, :, 1], 1, axis=1)
    X[0, :] = A[0, :, 0]
    X[-1, :] = -A[-2, :, 0]
    Y[:, 0] = A[:, 0, 1]
    Y[:, -1] = -A[:, -2, 1]
    return X + Y
    
"""
def perform_one_iteration(start, g, l, tau):
    expr = grad(div(start)-g/l)
    return (start + tau*expr)/(1+tau*np.abs(expr))
"""

def perform_one_iteration(start, g, l, tau):
    expr = grad(div(start)-g/l)
    return (start + tau*expr)/(1+tau*np.sqrt(expr[:,:, 0]**2 + expr[:,:, 1]**2))[:,:,np.newaxis]


def project(g, l, tau):
    current = np.zeros((*g.shape, 2))
    previous = np.full((*g.shape, 2), np.inf)
    while np.max(np.abs(current-previous))>0.01:        #0.01
        previous = current
        current = perform_one_iteration(current, g, l, tau)

    return div(current)*l


def denoise(g, sigma, tau=1/8, K = 35):
    N = g.shape[1]
    l = 1
    for _ in range(K):
        v = project(g, l, tau)
        f = np.linalg.norm(v)
        l = (N*sigma/f)*l
        print("Lambda:", l)
    return g-v

if __name__ == "__main__":
    sigma = 25
    test_image = 255*plt.imread("ethz_main.png")[:,:, 0]
    g = test_image + np.random.normal(size = test_image.shape, loc= 0, scale=sigma)
    cleaned = denoise(g, sigma)

    fig, ax = plt.subplots(ncols = 2)
    ax[0].axis('off')
    ax[1].axis('off')
    ax[0].imshow(g, cmap="gray")
    ax[0].set_title("Noisy", y=-0.12)
    ax[1].imshow(cleaned, cmap="gray")
    ax[1].set_title("Denoised", y=-0.12)
    fig.tight_layout()
    #plt.savefig("denoised_img.png", dpi=1200)
    plt.savefig("denoised_img.pdf", dpi=1200, bbox_inches='tight')
    plt.show()


