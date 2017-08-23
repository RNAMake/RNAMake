#coding=UTF-8
"""This file holds several functions for envisaging different data"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import auxi as a
import pdb
#import auxi as a

def envis6d(phi, hold_on=False):  # phi(x1,x2,x3,alpha,beta,gamma)
    """
    this will envisage the configuration of the SE(3) group
    """
    phi = np.abs(phi)
    threshold = (phi.max() - phi.min()) * 0.1 + phi.min()
    grid_phi = np.mgrid[0:phi.shape[0], 0:phi.shape[1], 0:phi.shape[2], 0:phi.shape[3], 0:np.pi:phi.shape[4] * 1j,
               0:phi.shape[5]]
    grid_phi = grid_phi.astype('float')
    grid_phi[3,:] *= 2*np.pi/phi.shape[3]
    grid_phi[5,:] *= 2*np.pi/phi.shape[5]
    grid_phi_list = grid_phi.reshape([6,-1])
    phi_list = phi.reshape([-1])
    effect_points = np.logical_and(phi_list >= threshold,phi_list > 0)

    effect_phi_list = phi_list[effect_points]
    effect_grid = grid_phi_list[:,effect_points]

    effect_grid_translation = effect_grid[:3,:]
    effect_grid_rotation = effect_grid[3:,:]
    effect_grid_rotation_representation = \
        lambda r: effect_grid_translation + \
                  r * np.array([
                      np.sin(effect_grid_rotation[1, :]) * np.sin(effect_grid_rotation[0, :]),
                      -np.sin(effect_grid_rotation[1, :]) * np.cos(effect_grid_rotation[0, :]),
                      np.cos(effect_grid_rotation[1, :])
                  ])
    effect_grid_spin_representation = \
        effect_grid_rotation_representation(0.5) + \
        0.05 * np.array([
            np.cos(effect_grid_rotation[2, :]) * np.cos(effect_grid_rotation[0, :]) -
            np.sin(effect_grid_rotation[2, :]) * np.cos(effect_grid_rotation[1, :]) * np.sin(
                effect_grid_rotation[0, :]),
            np.cos(effect_grid_rotation[2, :]) * np.sin(effect_grid_rotation[0, :]) +
            np.sin(effect_grid_rotation[2, :]) * np.cos(effect_grid_rotation[1, :]) * np.cos(
                effect_grid_rotation[0, :]),
            np.sin(effect_grid_rotation[2, :]) * np.sin(effect_grid_rotation[1, :])
        ])
    fig1 = plt.figure()
    fig1.canvas.set_window_title('SE(3) Illustration')
    ax = fig1.add_subplot(111,projection='3d')
    # ax.set_aspect('equal')
    ax.axis(xmin=-1,xmax=phi.shape[0],ymin=-1,ymax=phi.shape[1])
    ax.set_zlim(-1,phi.shape[2])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.scatter(effect_grid_translation[0], effect_grid_translation[1], effect_grid_translation[2], s=100,
               depthshade=False)
    # for trans in np.mgrid[0:0.6:3j]:
    #     ax.scatter(effect_grid_translation[0]+trans, effect_grid_translation[1], effect_grid_translation[2],c='c',marker='+')
    #     ax.scatter(effect_grid_translation[0], effect_grid_translation[1]+trans, effect_grid_translation[2], c='m',marker='+')
    #     ax.scatter(effect_grid_translation[0], effect_grid_translation[1], effect_grid_translation[2]+trans, c='y',marker='+')
    # ax.scatter(effect_grid_translation[0]+0.15, effect_grid_translation[1], effect_grid_translation[2]+0.5, c='r',marker='s')
    # for r in np.mgrid[0:0.4:5j]:
    #     ax.scatter(effect_grid_rotation_representation(r)[0],effect_grid_rotation_representation(r)[1],effect_grid_rotation_representation(r)[2],marker='.',c='r')
    ax.scatter(effect_grid_rotation_representation(0.5)[0],
               effect_grid_rotation_representation(0.5)[1],
               effect_grid_rotation_representation(0.5)[2], marker='*', c='r', depthshade=False)
    ax.scatter(effect_grid_spin_representation[0], effect_grid_spin_representation[1],
               effect_grid_spin_representation[2], marker='s', c=effect_phi_list, depthshade=False)
    if not hold_on:
        plt.show()


def envis3d(phi, hold_on=False):  #phi(x1,x2,x3)
    phi = np.abs(phi.copy())
    grid_cart = np.mgrid[0:phi.shape[0],0:phi.shape[1],0:phi.shape[2]]
    effect_points_list = (phi>0).reshape([-1])
    effect_phi = phi.reshape([-1])[effect_points_list]
    grid_cart_list = grid_cart.reshape([3,-1])
    effect_grid_list = grid_cart_list[:,effect_points_list]

    fig1 = plt.figure()
    fig1.canvas.set_window_title('Cartesian Plot')
    ax = fig1.add_subplot(111, projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.scatter(effect_grid_list[0], effect_grid_list[1], effect_grid_list[2], c=effect_phi, cmap='winter',
               depthshade=False)
    if not hold_on:
        plt.show()


def envis3dsph(phi_in, hold_on=False, y_neg=False, component='all', pos=True, tp=False):  # phi(r,phi,theta)
    if component == 'mag':
        phi = abs(phi_in)
    elif component == 'rel':
        phi = np.real(phi_in)
    elif component == 'img':
        phi = np.imag(phi_in)
    elif component == 'all':
        envis3dsph(phi_in, hold_on=True, y_neg=y_neg, component='mag', pos=pos, tp=tp)
        envis3dsph(phi_in, hold_on=True, y_neg=y_neg, component='rel', pos=pos, tp=tp)
        envis3dsph(phi_in, hold_on, y_neg, component='img', pos=pos, tp=tp)
        return
    else:
        return
    if tp:
        phi = phi.transpose([0, 2, 1])
    if y_neg:
        grid_sph_struc = np.mgrid[0:phi.shape[0], 0:2 * np.pi * (1.0 - 1.0 / phi.shape[1]):phi.shape[1] * 1j,
                         0:np.pi:phi.shape[2] * 1j]
        grid_sph = grid_sph_struc.reshape([3, -1]).astype('float')
    else:
        grid_sph = np.mgrid[0:phi.shape[0], 0:phi.shape[1], 0:phi.shape[2]].reshape([3, -1]).astype('float')
        grid_sph[1, :] *= 2 * np.pi / phi.shape[1]
        grid_sph[2, :] *= np.pi / phi.shape[2]
    if pos:
        effect_points = phi.reshape([-1]) > 0
        effect_points_list = phi.reshape([-1])[effect_points]
        grid_sph = grid_sph[:, effect_points]
    else:
        effect_points_list = phi.reshape([-1])
    grid_cart = np.array([grid_sph[0]*np.cos(grid_sph[1])*np.sin(grid_sph[2]),grid_sph[0]*np.sin(grid_sph[1])*np.sin(grid_sph[2]),grid_sph[0]*np.cos(grid_sph[2])])
    fig1 = plt.figure()
    fig1.canvas.set_window_title('Spherical Plot')
    if component == 'rel':
        fig1.canvas.set_window_title('Spherical Plot REAL PART')
    if component == 'img':
        fig1.canvas.set_window_title('Spherical Plot IMAG PART')
    ax = fig1.add_subplot(111, projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.scatter(grid_cart[0], grid_cart[1], grid_cart[2], c=effect_points_list, depthshade=False, cmap='winter')
    if not hold_on:
        plt.show()
