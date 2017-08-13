# coding=UTF-8

import numpy as np
import pdb

import auxi as a


def spft3_out(phi, flnm):  # phi(x1,x2,x3,alpha,beta,gamma)
    assert np.all(
        np.array(phi.shape[0:3]) == np.array(a.org) * 2 + 1)  # Should be odd * odd * odd and the origin at the center,
    # don't know what will happen otherwise.
    shape_out = np.array(phi.shape)
    phi_out = np.concatenate([shape_out, phi.reshape(-1)])
    np.savetxt(flnm, phi_out, fmt='%s+%sj')


def spft3_in(flnm):
    ifl = np.char.replace(np.loadtxt(flnm, delimiter=',', dtype=np.str), 'i', 'j').astype('complex')
    ifl = ifl.reshape(-1)
    ishape = ifl[0:6]
    idata = ifl[6:]
    phi = idata.reshape(np.real(ishape))
    a.rang_theta = phi.shape[2]
    return phi  # (r, phi, theta, alpha, beta, gamma)


def ispft3_out(phi, flnm):  # phi(a_p, theta, phi, lp, mp, s, n)
    print phi.shape
    assert phi.shape[1] % 2 == 1 and phi.shape[2] % 4 == 0
    shape_out = np.array(phi.shape)
    phi_out = np.concatenate([shape_out, phi.reshape(-1)])
    np.savetxt(flnm, phi_out, fmt='%.5e%+.5ej')


def ispft3_in(flnm):
    ifl = np.char.replace(np.loadtxt(flnm, delimiter=',', dtype=np.str), 'i', 'j').astype('complex')
    ifl = ifl.reshape(-1)
    ishape = ifl[0:7]
    idata = ifl[7:]
    phi = idata.reshape(np.real(ishape))
    return phi  # phi(x1,x2,x3,lp,mp,s,n)
