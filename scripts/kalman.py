from pyrat.diff import kalman as ka
import pyrat.diff.intfun as intfun
import numpy as np
import matplotlib.pyplot as plt
import itertools as it
import scipy as sp
import numpy as np
import pyrat.fileutils.gpri_files as gpf
import os
import matplotlib.colors as cols



def F_model(dt):
    #F matrix for variable dt
    F_m = np.array([[1, dt], [0, 1]])
    return F_m




def H_stack(f_fun, H_model, itab, t_vector):
    F_aug = []
    A = intfun.itab_to_incidence_matrix(itab)
    print(itab)
    print(A)
    F = np.eye(2)
    t_start  = t_vector[0]
    for idx_t, t in enumerate(t_vector[::]):
        t_end = t_vector[idx_t]
        dt = t_end - t_start
        F_model = f_fun(dt)
        F = np.dot(F_model,F)
        F_aug.append(np.dot(H_model, F))
        t_start = t_end
    F_aug = np.vstack(F_aug)
    pi = np.dot(A, F_aug)
    return pi

def plot_state_and_variance(ax, t, x_model, x_filter, P, **kawargs):
    ax.plot(t, x_filter)
    ax.plot(t, x_model)
    ax.fill_between(t, x_filter - P**0.5, y2=x_filter + P**0.5, facecolor='none')


n_stack  = 10
win = 3
stride = 2
step = 2
#Compute itab
itab = intfun.itab(n_stack, win, step, stride, 0)

# #Displacement model
lam = 3e8/17.9e9
H_d = np.array([[4 * np.pi/lam,0]])#Extracts the displacement


def kalman(input, output, threads, config, params, wildcards):
    #Number of stack epochs
    n_epochs = len(input.slc_par_names) // n_stack
    print(n_epochs)

    #All stack epochs
    for idx_epoch in range(1,n_epochs):
        # Create vector of times
        acquisition_times = []
        for idx_slc in range(idx_epoch, idx_epoch + n_stack):
            slc_par = gpf.par_to_dict(input.slc_par_names[idx_slc])
            acquisition_times.append(slc_par.start_time)
        H = H_stack(F_model, H_d, itab, acquisition_times)
        plt.plot(H[:,1])
        plt.show()
        print(H)





# #Covariance for the filter
# R_m = correlation_to_covariance(stack_correlation_matrix(itab, dt, sigma=sigma, gamma=0.8), np.diag((1e2,)*len(itab)))
# print(R_m)
# print(sigma)
# plt.imshow(R_m)
# plt.show()
# Q_f = np.eye(2) * 1e-2
# #Covariances for the model used to generate observations
# Q_m = np.zeros(F_m.shape)
# R_f = correlation_to_covariance(stack_correlation_matrix(itab, dt, sigma=sigma, gamma=0.6), np.diag((5e-1,)*len(itab)))
# #Starting vector for model
# x0_m = np.array([0,v])
# model =  ka.LinearSystem(F=F_m, H=H, Q=Q_m, R=R_m,x0=x0_m)
# #Construct kalman filter
# filter = ka.KalmanFilter(nstates=2, noutpus=H.shape[0],F=F_m,H=H, Q=Q_f, R=R_f, x0=[0,0.002])
#
#
# n_epochs = 50
# x_model = []
# x_filter = []
# z = []
# P = []
# #Simulate displacements
# for i in range(n_epochs):
#     x_model.append(model.state_transition())
#     z.append(model.output())
#     filter.predict()
#     filter.update(model.z)
#     x_filter.append(filter.x)
#     P.append(filter._P)
#
# x_filter = np.array(x_filter)
# x_model = np.array(x_model)
# z = np.array(z)
# P = np.array(P)
# t = np.arange(n_epochs) * dt
#
# f, ax_grid = plt.subplots(3,1, sharex=True)
# disp_ax = ax_grid[0]
# v_ax = ax_grid[1]
# plot_state_and_variance(disp_ax, t, x_model[:,0], x_filter[:,0], P[:,0,0])
# disp_ax.yaxis.set_label_text(r'Displacement [m]')
# disp_ax.xaxis.set_label_text(r'Time [s]')
# plot_state_and_variance(v_ax, t, x_model[:,1], x_filter[:,1], P[:,1,1])
# v_ax.yaxis.set_label_text(r'Velocity [m/s]')
# v_ax.set_ylim([-1e-2,1e-2])
# map = plt.get_cmap('RdBu')
# cNorm = cols.Normalize(vmin=0, vmax=n_slc)
# scalarMap = plt.cm.ScalarMappable(norm=cNorm, cmap=map)
# # for i in range(z.shape[0]):
# #     ax_grid[-1].plot(t, z[:, i], color=scalarMap.to_rgba(i), marker='o')
# plt.show()



kalman(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
      snakemake.wildcards)


