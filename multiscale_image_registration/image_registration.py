import numpy as np
from scipy.integrate import simps, dblquad
from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from itertools import product


def rk4(x0, N, f, hook = lambda x: x):
        xs = np.zeros((N, *x0.shape))
        ts, h = np.linspace(0, 1, N, retstep = True)
        xs[0] = np.copy(x0)
        for i in range(1, N):
            k1 = f(ts[i-1], xs[i-1])
            k2 = f(ts[i-1] + h/2, xs[i-1] + k1*h/2)
            k3 = f(ts[i-1] + h/2, xs[i-1] + k2*h/2)
            k4 = f(ts[i-1] + h, xs[i-1] + k3*h)
            xs[i] = hook(xs[i-1] + (k1 + 2*k2 + 2*k3+ k4)*h/6)

        return ts, xs


shape_functions ={0: lambda x: (2*np.abs(x)**3 - 3*np.abs(x)**2 + 1)*(np.abs(x) < 1),
                  1: lambda x: np.sign(x)*(np.abs(x)**3 - 2*np.abs(x)**2 + np.abs(x))*(np.abs(x) < 1)}

shape_functions_der ={0: lambda x: np.sign(x)*(6*np.abs(x)**2 - 6*np.abs(x)),
                  1: lambda x: 3*np.abs(x)**2 - 4*np.abs(x) + 1,}


class FiniteElementFunction:
    def __init__(self, n_t, n_x):
        self.n_t = n_t
        self.n_x = n_x
        self.parameters = np.zeros((n_t, n_x, 2))
        self.x_supports, self.dx = np.linspace(-0.000001, 1.000001, n_x, retstep=True)
        self.t_supports, self.dt = np.linspace(-0.000001, 1.000001, n_t, retstep=True)

        self.l2_dict = {}

        for i, j in product(range(0, 2), range(0, 2)):
            shape_i = lambda t, x: t/self.dt*shape_functions[i](x/self.dx)
            shape_j= lambda t, x: t/self.dt*shape_functions[j](x/self.dx)
            shape_i_translate_t = lambda t, x: (self.dt - t)/self.dt*shape_functions[i](x/self.dx)
            shape_i_translate_x = lambda t, x: t/self.dt*shape_functions[i]((self.dx - x)/self.dx)
            shape_i_translate_tx = lambda t, x: (self.dt - t)/self.dt*shape_functions[i]((self.dx - x)/self.dx)
            
            self.l2_dict[i, j, 0, 0] = 2*dblquad(lambda t, x: shape_i(t, x)*shape_j(t, x), -self.dx, self.dx, 0, self.dt)[0]
            self.l2_dict[i, j, 1, 0] = dblquad(lambda t, x: shape_i_translate_t(t, x)*shape_j(t, x), -self.dx, self.dx, 0, self.dt)[0]
            self.l2_dict[i, j, 1, 1] = dblquad(lambda t, x: shape_i_translate_tx(t, x)*shape_j(t, x), -self.dx, self.dx, 0, self.dt)[0]
            self.l2_dict[i, j, 0, 1] = 2*dblquad(lambda t, x: shape_i_translate_x(t, x)*shape_j(t, x), -self.dx, self.dx, 0, self.dt)[0]

        self.flat_to_ind = {}
        self.ind_to_flat = {}
        counter = 0
        
        for ind in np.ndindex((n_t, n_x, 2)):
            if ind[1] != 0 and ind[1] != n_x - 1:
                self.flat_to_ind[counter] = ind
                self.ind_to_flat[ind] = counter
                counter += 1

    def evaluate(self, t, x):
        x = np.clip(x, 0, 1)
        t = np.clip(t, 0, 1)
        time_bin = np.digitize(t, self.t_supports)
        space_bin = np.digitize(x, self.x_supports)

        if isinstance(x, np.ndarray):
            result = np.zeros_like(x)
        else:
            result = 0
        

        result += (t - self.t_supports[time_bin-1])/self.dt * self.parameters[time_bin, space_bin-1, 0]*shape_functions[0]((x - self.x_supports[space_bin - 1])/self.dx)
        result += (t - self.t_supports[time_bin-1])/self.dt * self.parameters[time_bin, space_bin, 0]*shape_functions[0]((x - self.x_supports[space_bin])/self.dx)

        result += (t - self.t_supports[time_bin-1])/self.dt * self.parameters[time_bin, space_bin-1, 1]*shape_functions[1]((x - self.x_supports[space_bin - 1])/self.dx)
        result += (t - self.t_supports[time_bin-1])/self.dt * self.parameters[time_bin, space_bin, 1]*shape_functions[1]((x - self.x_supports[space_bin])/self.dx)
        
        result += (self.t_supports[time_bin] - t)/self.dt * self.parameters[time_bin - 1, space_bin-1, 0]*shape_functions[0]((x - self.x_supports[space_bin - 1])/self.dx)
        result += (self.t_supports[time_bin] - t)/self.dt * self.parameters[time_bin - 1, space_bin, 0]*shape_functions[0]((x - self.x_supports[space_bin])/self.dx)

        result += (self.t_supports[time_bin] - t)/self.dt * self.parameters[time_bin - 1, space_bin-1, 1]*shape_functions[1]((x - self.x_supports[space_bin - 1])/self.dx)
        result += (self.t_supports[time_bin] - t)/self.dt * self.parameters[time_bin - 1, space_bin, 1]*shape_functions[1]((x - self.x_supports[space_bin])/self.dx)
        return result

    def derive(self, t, x):
        time_bin = np.digitize(t, self.t_supports)
        space_bin = np.digitize(x, self.x_supports)

        if isinstance(x, np.ndarray):
            result = np.zeros_like(x)
        else:
            result = 0
        

        result += (t - self.t_supports[time_bin-1])/self.dt * self.parameters[time_bin, space_bin-1, 0]*shape_functions_der[0]((x - self.x_supports[space_bin - 1])/self.dx)/self.dx
        result += (t - self.t_supports[time_bin-1])/self.dt * self.parameters[time_bin, space_bin, 0]*shape_functions_der[0]((x - self.x_supports[space_bin])/self.dx)/self.dx

        result += (t - self.t_supports[time_bin-1])/self.dt * self.parameters[time_bin, space_bin-1, 1]*shape_functions_der[1]((x - self.x_supports[space_bin - 1])/self.dx)/self.dx
        result += (t - self.t_supports[time_bin-1])/self.dt * self.parameters[time_bin, space_bin, 1]*shape_functions_der[1]((x - self.x_supports[space_bin])/self.dx)/self.dx
        
        result += (self.t_supports[time_bin] - t)/self.dt * self.parameters[time_bin - 1, space_bin-1, 0]*shape_functions_der[0]((x - self.x_supports[space_bin - 1])/self.dx)/self.dx
        result += (self.t_supports[time_bin] - t)/self.dt * self.parameters[time_bin - 1, space_bin, 0]*shape_functions_der[0]((x - self.x_supports[space_bin])/self.dx)/self.dx

        result += (self.t_supports[time_bin] - t)/self.dt * self.parameters[time_bin - 1, space_bin-1, 1]*shape_functions_der[1]((x - self.x_supports[space_bin - 1])/self.dx)/self.dx
        result += (self.t_supports[time_bin] - t)/self.dt * self.parameters[time_bin - 1, space_bin, 1]*shape_functions_der[1]((x - self.x_supports[space_bin])/self.dx)/self.dx
        return result

    def l2_product(self, ind_i, ind_j):
        i_t, i_x, i_f = ind_i
        j_t, j_x, j_f = ind_j
        
        if abs(i_t - j_t) > 1 or abs(i_x - j_x)>1:
            return 0

        result =  self.l2_dict[(i_f, j_f, abs(i_t - j_t), abs(i_x - j_x))]

        if i_t == j_t and (i_t == 0 or i_t == self.n_t - 1):
            return 0.5*result
        return result

    def compute_stiffness_matrix(self):
        columns = []
        rows = []
        data = []
        for ind_a, ind_b in product(self.ind_to_flat, self.ind_to_flat):
            entry = self.l2_product(ind_a, ind_b)
            if entry != 0:
                columns.append(self.ind_to_flat[ind_a])
                rows.append(self.ind_to_flat[ind_b])
                data.append(entry)

        self.stiffness_matrix = csc_matrix((data, (rows, columns)), shape = (len(self.ind_to_flat), len(self.ind_to_flat)))
        




class VelocityField(FiniteElementFunction):
    def __init__(self, n_t, n_x, N_t, N_x, I0):
        super().__init__(n_t, n_x)
        self.I0 = I0
        self.flow_N_t = N_t
        self.flow_N_x = N_x

        self.flow_x_supports, self.flow_dx = np.linspace(0, 1, self.flow_N_x, retstep=True)
        self.flow_t_supports, self.flow_dt = np.linspace(0, 1, self.flow_N_t, retstep=True)

        self.grad_I0 = np.zeros_like(self.I0)
        self.grad_I0[1:-1] = (self.I0[2:] - self.I0[:-2])/(2*(len(self.I0) - 1))
        self.grad_I0[0] = (self.I0[1] - self.I0[0])/(len(self.I0) - 1)
        self.grad_I0[-1] = (self.I0[-2] - self.I0[-1])/(len(self.I0) - 1)
        self.image_supports = np.linspace(0, 1, len(self.I0))

    def compute_flow(self):
        self.flow = rk4(self.flow_x_supports, self.flow_N_t, self.evaluate, hook = lambda x: np.clip(x, 0, 1))[1]
        self.transformed_grad_I0 = np.interp(self.flow[-1], image_supports, self.grad_I0)
        
    def compute_mu(self):
        self.mu = rk4(np.ones_like(self.flow_x_supports), self.flow_N_t, lambda t, x: self.derive(t, self.flow[np.rint(t/self.dt).astype(int)])*x)[1]
        
    def frechet_derivative(self, df):
        expanded_ts = np.repeat(np.expand_dims(self.flow_t_supports, 1), self.flow_N_x, axis = 1)
        
        return self.transformed_grad_I0 * self.mu[-1]*simps(1/self.mu * df.evaluate(expanded_ts, self.flow), axis = 0)

    def frechet_adjoint(self, I):
        print("Adjoint")
        rhs = np.zeros(len(self.ind_to_flat))
        
        basis_element = FiniteElementFunction(self.n_t, self.n_x)

        for ind in self.ind_to_flat:
            basis_element.parameters[ind] = 1
            img = self.frechet_derivative(basis_element)
            basis_element.parameters[ind] = 0
            rhs[self.ind_to_flat[ind]] = np.sum(img*I)/len(I)
        
        result = spsolve(self.stiffness_matrix, rhs)
        reshaped = np.zeros_like(self.parameters)
        for ind in self.ind_to_flat:
            reshaped[ind] = result[self.ind_to_flat[ind]]

        return reshaped



def warp(x):
   return x - 0.125*np.cos(9*(2*x-1))*np.exp(-1/(1-(2*x-1)**2))


def ista_optimize(I0, I1, l, step_size, n_t, n_x, N_t):
    N_x = len(I0)
    image_supports = np.linspace(0, 1, len(I0))
    vel_field = VelocityField(n_t, n_x, N_t, N_x, I0)
    vel_field.compute_stiffness_matrix()
    for _ in range(80):
        vel_field.compute_flow()
        vel_field.compute_mu()
        current_image = np.interp(vel_field.flow[-1], image_supports, I0)
        vel_field.parameters = shrinkage_operator(vel_field.parameters - 2*step_size*vel_field.frechet_adjoint(current_image - I1),
                           l*step_size)

    vel_field.compute_flow()
    return vel_field


def shrinkage_operator(x, alpha):
    return np.maximum((np.abs(x) - alpha), 0)*np.sign(x)


def multiscale(I0, I1, lambdas, step_size, n_t, n_x, N_t):
    current_image = np.copy(I0)
    image_supports = np.linspace(0, 1, len(I0))
    flows = []
    for l in lambdas:
        vel_field = ista_optimize(current_image, I1, l, step_size, n_t, n_x, N_t)
        flows.append(np.copy(vel_field.flow))
        current_image = np.interp(vel_field.flow[-1], image_supports, current_image)
    return flows

def compose_flows(flows):
    compositions = [np.copy(flows[0][-1])]
    for next_flow in flows[1:]:
        compositions.append(np.interp(next_flow[-1], next_flow[0], compositions[-1]))
    return compositions

def animate(I0, I1, compositions, file_name = None):
    fig,ax = plt.subplots(1)
    line1, = ax.plot(compositions[0], I1)
    line2, = ax.plot(compositions[0], np.interp(compositions[0], compositions[0], I0))
    def animation_function(i):
        line2.set_ydata(np.interp(compositions[i], compositions[0], I0))
        return line1, line2
    def init_func():
        line2.set_ydata([np.nan]*len(compositions[0]))
        return line1, line2

    ani = animation.FuncAnimation(fig, animation_function, init_func = init_func, frames = len(compositions), interval=5000/len(compositions), blit=True)
    if file_name is not None:
        ani.save(file_name)
    plt.show()

def compose_full(flows):
    compositions = np.zeros((len(flows)*(flows[0].shape[0] - 1) + 1, flows[0].shape[1]))
    compositions[0:flows[0].shape[0]] = flows[0]
    current_anchor = flows[0].shape[0] - 1
    i = flows[0].shape[0]
    for flow in flows[1:]:
        for step in flow[1:]:
            compositions[i] = np.interp(step, flow[0], compositions[current_anchor])
            i+=1
        current_anchor = i - 1
    return compositions

if __name__ == "__main__":
    n_t = 50
    n_x = 5             #5
    N_t = 100           #100
    N_x = 100           #100
    start = 0
    image_supports = np.linspace(0, 1, N_x)
    image = np.sin(8*np.pi*image_supports**2)
    transformed_image = np.interp(warp(image_supports), image_supports, image)
    single = False

    if single:
        vel_field = ista_optimize(transformed_image, image, 0., 0.005, n_t, n_x, N_t)
        vel_field.compute_flow()
        plt.plot(vel_field.flow[-1])
        plt.plot(warp(image_supports))
        plt.show()

        print(f"Original L2 distance {np.sum((image - transformed_image)**2)}, now {np.sum((image - np.interp(vel_field.flow[-1], image_supports, transformed_image))**2)}")
        plt.plot(np.interp(vel_field.flow[-1], image_supports, transformed_image))
        plt.show()
    else:
        flows = multiscale(transformed_image, image, (3.**(-np.arange(start, start+5))), 0.01, n_t, n_x, N_t)                    #0.01
        compositions = compose_flows(flows)
        fig, ax = plt.subplots(nrows = len(compositions) + 1)
        ax[0].plot(image_supports, image, label="$I_1$")
        ax[0].plot(image_supports, transformed_image, label = "$I_0$")
        ax[0].xaxis.set_ticklabels([])
        ax[0].yaxis.set_ticklabels([])
        ax[0].legend(loc='upper right')
        for i, composition in enumerate(compositions):
                ax[i+1].plot(image_supports, image, label="$I_1$")
                comp_string = f"\\circ g_0 \\circ \\dots \\circ g_{i}" if i>2 else "\\circ " + "\\circ ".join([f"g_{j}" for j in range(0, i+1)])
                ax[i+1].plot(image_supports, np.interp(composition, image_supports, transformed_image), label = f"$I_0 {comp_string}$")
                ax[i+1].legend(loc='upper right')
                ax[i+1].xaxis.set_ticklabels([])
                ax[i+1].yaxis.set_ticklabels([])
                ax[i+1].set_ylabel(f"$\\lambda = 3^{i+start}$", rotation=0, labelpad=20)
        fig.tight_layout()
        plt.savefig("multiscale_decomposition.pdf", dpi=1200, bbox_inches='tight')
        plt.show()

        full_compositions = compose_full(flows)
        animate(transformed_image, image, full_compositions, "morphing.mp4")
