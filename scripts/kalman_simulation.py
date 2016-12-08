from pyrat.diff import kalman as ka
import pyrat.diff.intfun as intfun
import numpy as np
import matplotlib.pyplot as plt
import itertools as it

dt = 5*60
v = 2/(24*60*60)
a = 0.0
n_slc = 8

#Standard deviation of atmosphere
lam = 3e8/17.9
atm_sigma = 1e4/(24*60*60)
sigma = atm_sigma * 4 * np.pi/lam
print(sigma)

#Displacement model
F_m = np.array([[1,dt],[0,1]])
H_d = np.array([[4 * np.pi/lam,0]])#Extracts the displacement

def H_stack(F_model, H_model, itab):
    H = []
    for master_idx, slave_idx, *rest in itab:
        H.append(np.dot(H_model, (np.linalg.matrix_power(F_model,master_idx) - np.linalg.matrix_power(F_model,slave_idx))))
    H = np.vstack(H)
    return H


def brownian_variance(gamma, sigma, dt):
    return gamma * np.exp(-dt/(2/(sigma**2)))

def interferogram_coherence(gamma, sigma, dt):
    return brownian_variance(gamma, sigma, dt) *np.exp(1j * (dt)*disp

#using roccas model
def stack_correlation_matrix(itab, dt, gamma=0.8, sigma=1e-3):
    C = np.eye(len(itab))
    for idx_stack, (master_idx, slave_idx, *rest) in enumerate(itab):
        for idx_stack_1, (master_idx_1, slave_idx_1, *rest) in enumerate(itab):
            dt1 = np.abs(master_idx - master_idx_1) * dt
            dt2 = np.abs(slave_idx - slave_idx_1) * dt
            c1 = interferogram_covariance(gamma, sigma, dt1)
            c2 = interferogram_covariance(gamma, sigma, dt2)
            print(master_idx, master_idx_1)
            C[idx_stack, idx_stack_1] =  (c1 * c2)
    return C


def correlation_to_covariance(C, p)

def plot_state_and_variance(ax, t, x_model, x_filter, P, **kawargs):
    ax.plot(t, x_filter)
    ax.plot(t, x_model)
    ax.fill_between(t, x_filter - P**0.5, y2=x_filter + P**0.5, facecolor='none')


#Compute itab
itab = intfun.itab(n_slc, n_slc, 1, 1, 0)
#Compute H matrix from Itab
H = H_stack(F_m, H_d, itab)
#Covariance for the filter
R_m = stack_correlation_matrix(itab, dt, sigma=sigma)
print(sigma)
plt.imshow(R_m)
plt.show()
Q_f = np.eye(2) * 1e-2
#Covariances for the model used to generate observations
Q_m = np.zeros(F_m.shape)
R_f = stack_correlation_matrix(itab, dt, sigma=sigma)
#Starting vector for model
x0_m = np.array([0,v])
model =  ka.LinearSystem(F=F_m, H=H, Q=Q_m, R=R_m,x0=x0_m)
#Construct kalman filter
filter = ka.KalmanFilter(nstates=2, noutpus=H.shape[0],F=F_m,H=H, Q=Q_f, R=R_f, x0=[0,0.002])


n_epochs = 50
x_model = []
x_filter = []
z = []
P = []
#Simulate displacements
for i in range(n_epochs):
    x_model.append(model.state_transition())
    z.append(model.output())
    filter.predict()
    filter.update(model.z)
    x_filter.append(filter.x)
    P.append(filter._P)

x_filter = np.array(x_filter)
x_model = np.array(x_model)
z = np.array(z)
P = np.array(P)
t = np.arange(n_epochs) * dt

import matplotlib.colors as cols
f, ax_grid = plt.subplots(3,1, sharex=True)
disp_ax = ax_grid[0]
v_ax = ax_grid[1]
plot_state_and_variance(disp_ax, t, x_model[:,0], x_filter[:,0], P[:,0,0])
disp_ax.yaxis.set_label_text(r'Displacement [m]')
disp_ax.xaxis.set_label_text(r'Time [s]')
plot_state_and_variance(v_ax, t, x_model[:,1], x_filter[:,1], P[:,1,1])
v_ax.yaxis.set_label_text(r'Velocity [m/s]')
v_ax.set_ylim([-5,5])
map = plt.get_cmap('RdBu')
cNorm = cols.Normalize(vmin=0, vmax=n_slc)
scalarMap = plt.cm.ScalarMappable(norm=cNorm, cmap=map)
# for i in range(z.shape[0]):
#     ax_grid[-1].plot(t, z[:, i], color=scalarMap.to_rgba(i), marker='o')
plt.show()


