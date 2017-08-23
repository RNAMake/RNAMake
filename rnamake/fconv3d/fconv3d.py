# coding=UTF-8
"""Main routines for convolution"""

import pdb

import numpy as np
import numpy.fft as ft
import scipy.interpolate as ip
from numba import jit
import warnings
# import tricubic

warnings.filterwarnings("error", category=np.ComplexWarning)

import auxi as a
import basis as bs
import numexpr as ne


def c3d(phi1, phi2):  # phi1(x1,x2,x3,alpha,beta,gamma):

    f11 = e_fft3(phi1)
    del phi1
    f12 = e_fft3(phi2)  # f1(p1,p2,p3,alpha,beta,gamma)
    del phi2

    # f1_sph(a_p,phi,theta,alpha,beta,gamma)
    f11_sph = car2sph(f11)
    del f11
    f12_sph = car2sph(f12)
    del f12

    # f2(a_p,phi,theta,lp,mp,n)
    f21 = so3int(f11_sph)
    del f11_sph
    f22 = so3int(f12_sph)
    del f12_sph

    # f_hat(a_p,m,lp,mp,s,l)
    f_hat1 = usphint(f21)
    del f21
    f_hat2 = usphint(f22)
    del f22

    f_hat = mltp(f_hat1, f_hat2)
    del f_hat1, f_hat2

    # f = inv_int(f_hat)
    #
    g1 = usphinti(f_hat)  # g1(a_p,theta,phi,lp,mp,s,n)

    g2 = sph2cart(g1)  # g2(p1,p2,p3,lp,mp,s,n)

    g3 = e_ifft3(g2)  # g3(x1,x2,x3,lp,mp,s,n)

    f = so3inti(g3)

    return f


def e_fft3(phi):
    print 'This is e_fft3 ... \n'
    f = np.copy(phi)
    f = np.concatenate([f[a.org[0]:, :, :, :, :, :], f[:a.org[0], :, :, :, :, :]], axis=0)
    f = np.concatenate([f[:, a.org[1]:, :, :, :, :], f[:, :a.org[1], :, :, :, :]], axis=1)
    f = np.concatenate([f[:, :, a.org[2]:, :, :, :], f[:, :, :a.org[2], :, :, :]], axis=2)
    res = ft.ifftn(f, axes=(0, 1, 2))
    '''normalization'''
    # res *= phi.shape[0] * phi.shape[1] * phi.shape[2]
    # res /= (2 * np.pi) ** 3
    res_shifted = ft.fftshift(res, axes=(0, 1, 2))
    print 'Goodbye e_fft3. \n'
    return res_shifted


@jit
def car2sph(phi):
    # phi(p1,p2,p3,alpha,beta,gamma)
    #     input expectations:
    # p1,p2,p3 same, odd
    # alpha, gamma same,odd
    # beta
    print  'This is car2sph... \n'
    f = phi
    rang_p1, rang_p2, rang_p3, rang_al, rang_be, rang_ga = f.shape
    a.rang_m = 4*(rang_p1-1)+1
    assert rang_p1 == rang_p2 == rang_p3 and rang_p1 % 2 == 1 and \
           rang_al == rang_ga and rang_al % 2 == 1
    a.rang_p = [rang_p1, rang_p2, rang_p3]
    p_max = (rang_p1 - 1) / 2
    p_min = -p_max
    rang_r = 4 * np.ceil((rang_p1 - 1) * 1.732 / 4) + 1
    factor = (rang_p1 - 1.0) / (rang_r - 1.0)
    a.factor = factor
    rang_rp = 4 * np.ceil((rang_p1 - 1) * 1.732 / 4) + 1
    rang_phi = np.ceil(rang_rp * 1.414 / 1.732 / 2) * 4
    if a.rang_m > rang_phi:
        a.rang_m = rang_phi -1
    a.rang_m = int(a.rang_m)
    rang_theta = np.ceil(rang_rp / 2) * 2 + 1
    axis_ticks = (np.mgrid[p_min:p_max + 1], np.mgrid[p_min:p_max + 1],
                  np.mgrid[p_min:p_max + 1])
    sph_grid = np.mgrid[0:rang_r, 0:rang_phi, 0:np.pi:rang_theta * 1j].astype('float32')
    sph_grid[1, :] *= 2 * np.pi / rang_phi
    r_g, phi_g, theta_g = sph_grid[0], sph_grid[1], sph_grid[2]
    cart_grid = np.array([r_g * np.cos(phi_g) * np.sin(theta_g),
                          r_g * np.sin(phi_g) * np.sin(theta_g),
                          r_g * np.cos(theta_g)]).transpose([1, 2, 3, 0]) * factor
    res = np.zeros([rang_r, rang_phi, rang_theta, rang_al, rang_be, rang_ga]).astype('complex64')
    for alpha in range(rang_al):
        for beta in range(rang_be):
            for gamma in range(rang_ga):
                res[:, :, :, alpha, beta, gamma] = \
                    ip.interpn(axis_ticks, f[:, :, :, alpha, beta, gamma], cart_grid,
                               bounds_error=False, fill_value=0)
    # res[0,1:,1:,:,:,:] = 0
    # res[:,1:,0,:,:,:] = 0
    # res[:,1:,-1,:,:,:] = 0

    print 'Goodbye car2sph. \n'
    return res


@jit
def car2sph3(phi):
    # phi(p1,p2,p3,alpha,beta,gamma)
    #     input expectations:
    # p1,p2,p3 same, odd
    # alpha, gamma same,odd
    # beta
    print  'This is car2sph3... \n'
    f = phi
    rang_p1, rang_p2, rang_p3, rang_al, rang_be, rang_ga = f.shape
    assert rang_p1 == rang_p2 == rang_p3 and rang_p1 % 2 == 1 and \
           rang_al == rang_ga and rang_al % 2 == 1
    a.rang_p = [rang_p1, rang_p2, rang_p3]
    p_max = (rang_p1 - 1) / 2
    p_min = -p_max
    rang_r = int(2 * np.ceil((rang_p1 - 1) * 1.732 / 4) + 1)
    factor = (rang_p1 - 1.0) / (rang_r - 1.0)
    a.factor = factor
    rang_rp = int(2 * np.ceil((rang_p1 - 1) / 4 * 1.732) + 1)
    rang_phi = int(np.ceil(rang_rp * 1.414 / 1.732 / 2) * 4)
    rang_theta = int(np.ceil(rang_rp / 2) * 2 + 1)
    # axis_ticks = (np.mgrid[p_min:p_max + 1], np.mgrid[p_min:p_max + 1],
    # np.mgrid[p_min:p_max + 1])
    sph_grid = np.mgrid[0:rang_r, 0:rang_phi, 0:np.pi:rang_theta * 1j].astype('float64')
    sph_grid[1, :] *= 2 * np.pi / rang_phi
    r_g, phi_g, theta_g = sph_grid[0], sph_grid[1], sph_grid[2]
    cart_grid = np.array([r_g * np.cos(phi_g) * np.sin(theta_g) - p_min,
                          r_g * np.sin(phi_g) * np.sin(theta_g) - p_min,
                          r_g * np.cos(theta_g) - p_min]).transpose([1, 2, 3, 0]) * factor
    res = np.zeros([rang_r, rang_phi, rang_theta, rang_al, rang_be, rang_ga]).astype('complex64')
    for alpha in range(rang_al):
        print alpha
        for beta in range(rang_be):
            for gamma in range(rang_ga):
                # print alpha, beta, gamma
                f_tmp = f[:, :, :, alpha, beta, gamma]
                f_tmp_r = np.real(f_tmp)
                f_tmp_i = np.imag(f_tmp)
                itp_r = tricubic.tricubic(list(f_tmp_r), [rang_p1, rang_p2, rang_p3])
                itp_i = tricubic.tricubic(list(f_tmp_i), [rang_p1, rang_p2, rang_p3])
                for a_p in range(rang_r):
                    for fi in range(rang_phi):
                        for theta in range(rang_theta):
                            res[a_p, fi, theta, alpha, beta, gamma] = itp_r.ip([cart_grid[a_p, fi, theta, 0],
                                                                                cart_grid[a_p, fi, theta, 1],
                                                                                cart_grid[a_p, fi, theta, 2]]) + \
                                                                      itp_i.ip([cart_grid[a_p, fi, theta, 0],
                                                                                cart_grid[a_p, fi, theta, 1],
                                                                                cart_grid[a_p, fi, theta, 2]]) * 1j

    print 'Goodbye car2sph3. \n'
    return res


@jit
def so3int(phi):
    # phi(a_p,phi,theta,alpha,beta,gamma)
    #     Input shape expectations
    # a_p
    # phi doubly even
    # theta odd
    # alpha odd
    # beta
    # gamma odd
    # alpha gamma same
    print 'This is so3int... \n'
    f = phi
    rang_ap, rang_phi, rang_theta, rang_al, rang_be, rang_ga = phi.shape
    '''dimension check'''
    assert rang_phi % 4 == 0 and rang_theta % 2 == 1 \
           and rang_al % 2 == 1 and rang_ga == rang_al
    a.n_max = (rang_al - 1) / 2
    a.n_min = - a.n_max
    rang_n = a.n_max - a.n_min + 1
    rang_mp = rang_n
    assert rang_n == rang_ga
    a.mp_max = (rang_ga - 1) / 2
    a.mp_min = - a.mp_max
    rang_lp = a.n_max + 1
    a.rang_beta = rang_be
    '''FFT over alpha and gamma'''
    f1 = ft.ifft(f, axis=3)
    f2 = ft.ifft(f1, axis=5)
    # del f1
    # f2(a_p,phi,theta,n,beta,mp)

    '''calculation of basis function'''
    # TODO values of contiguous indices of basis function can be generated much faster by recurrence formulae,
    # TODO but for now this step is not the bottleneck so that is not used.
    p = np.mgrid[0:rang_lp, 0:rang_al, 0:rang_ga, 0:rang_be].astype('float32')
    # p1(lp,n,mp,beta)
    p[3, :] += 0.5
    p[3, :] = np.cos(p[3, :] * np.pi / rang_be)

    p[1, :] += a.n_min
    if np.any(p[1, :, -a.n_min, :, :] != 0):
        raise 'wrong calculation!'
    p[1, :] = np.concatenate([p[1, :, -a.n_min:, :, :],
                              p[1, :, :-a.n_min, :, :]], axis=1)

    p[2, :] += a.mp_min
    p[2, :] = np.concatenate([p[2, :, :, -a.mp_min:, :],
                              p[2, :, :, :-a.mp_min, :]], axis=2)
    p1 = bs.p_lmn(p[0, :], p[1, :], p[2, :], p[3, :]) * (-1) ** p[1, :] * (-1) ** p[2, :]
    # del p
    # p1 = bs.p_lmn(p[0, :], p[1, :], p[2, :], p[3, :])
    # p1(lp,n,mp,beta)

    '''measure or weight'''
    # TODO here a proper weight for quadrature rather than simply measure for integration is needed.
    sinbeta = np.sin((np.mgrid[0:rang_be].astype('float32') + 0.5) * np.pi / rang_be)[np.newaxis, np.newaxis,
              np.newaxis,
              np.newaxis, np.newaxis, np.newaxis,
              :]

    '''sum to answer'''
    p2 = p1[np.newaxis, np.newaxis, np.newaxis, :, :, :, :]
    # p2(a_p,phi,theta,lp,n,mp,beta)

    # f2(a_p,phi,theta,n,beta,mp)
    # f2a(a_p,phi,theta,lp,n,mp,beta)
    f2a = f2.transpose([0, 1, 2, 3, 5, 4])[:, :, :, np.newaxis, :, :, :]
    res = np.zeros([rang_ap, rang_phi, rang_theta, rang_lp, rang_n, rang_mp]).astype('complex64')
    for ap in range(rang_ap):
        for fi in range(rang_phi):
            for theta in range(rang_theta):
                for lp in range(rang_lp):
                    for n in range(rang_n):
                        res[ap, fi, theta, lp, n, :] = \
                            np.trapz(
                                p2[0, 0, 0, lp, n, :, :] * f2a[ap, fi, theta, 0, n, :, :] * sinbeta[0, 0, 0, 0, 0, :,
                                                                                            :], dx=(np.pi / rang_be),
                                axis=1) + \
                            p2[0, 0, 0, lp, n, :, 0] * f2a[ap, fi, theta, 0, n, :, 0] * sinbeta[0, 0, 0, 0, 0, :,
                                                                                        0] * np.pi / rang_be / 4 + \
                            p2[0, 0, 0, lp, n, :, -1] * f2a[ap, fi, theta, 0, n, :, -1] * sinbeta[0, 0, 0, 0, 0, :,
                                                                                          -1] * np.pi / rang_be / 4
    # res = np.sum(p2 * f2a * sinbeta, axis=6)
    # res(a_p,phi,theta,lp,n,mp)
    resa = res.transpose([0, 1, 2, 3, 5, 4])
    # resa(a_p,phi,theta,lp,mp,n)

    print 'Goodbye so3int. \n'
    return resa


@jit
def usphint(phi):
    # phi(a_p,phi,theta,lp,mp,n)
    #     input expectations:
    # a_p
    # phi doubly even
    # theta odd
    # lp (mp+1)/2
    # mp odd
    # n odd
    # mp n same
    print 'This is usphint... \n'
    f = phi  # f(a_p,phi,theta,lp,mp,n)
    rang_ap, rang_phi, rang_theta, rang_lp, rang_mp, rang_n = f.shape
    a.rang_theta = rang_theta
    a.rang_phi = rang_phi

    '''dimension check'''
    assert rang_phi % 4 == 0 and rang_theta % 2 == 1 \
           and rang_mp % 2 == 1 and rang_mp == rang_n == rang_lp * 2 - 1

    '''integral over phi'''
    rang_s = rang_mp
    a.s_min = a.mp_min
    a.s_max = a.mp_max
    rang_m = a.rang_m
    a.m_max = (rang_m-1) / 2
    a.m_min = -a.m_max
    rang_l = (a.m_max+1)
    nphi = np.mgrid[0:rang_phi, 0:rang_n].astype('complex64')
    nphi[1, :] += a.n_min
    nphi[1, :] = np.concatenate([nphi[1, :, -a.n_min:],
                                 nphi[1, :, :-a.n_min]], axis=1)
    nphi[0, :] *= 2 * np.pi / rang_phi
    eptnphi = np.exp(-1j * nphi[1, :] * nphi[0, :])
    f2e = f * eptnphi[np.newaxis, :, np.newaxis, np.newaxis, np.newaxis, :]
    f3m = ft.ifft(f2e, axis=1)  # f3(a_p,m,theta,lp,mp,n)
    # del f2e


    '''truncating frequency domain'''
    f3m = np.concatenate([f3m[:,:a.m_max+1,:,:,:,:],f3m[:,a.m_min:,:,:,:,:]],axis = 1)

    # f3 = ft.ifft(f, axis=1) * rang_phi  # f3(a_p,m-n,theta,lp,mp,n)
    # rang_mn = rang_phi - 1
    # m_n_max = rang_phi / 2 - 1
    # m_n_min = - m_n_max
    # rang_m = rang_mn + rang_n - 1
    # print 'rang_m', rang_m
    # print 'rang_mn', rang_mn
    # print 'rang_n', rang_n
    # a.m_max = m_n_max + a.n_max
    # a.m_min = m_n_min + a.n_min
    #
    # ''''conversation of subscript'''
    # f3m = np.zeros(
    #     [rang_ap, rang_mn + rang_n - 1, rang_theta, rang_lp, rang_mp, rang_n]) + 0j
    # for m_n in xrange(m_n_min, m_n_max + 1, 1):
    #     for n in range(a.n_min, a.n_max + 1, 1):
    #         f3m[:, m_n + n, :, :, :, n] += f3[:, m_n, :, :, :, n]
    # # f3m(a_p,m,theta,lp,mp,n)
    # del f3

    '''calculation of basis function'''
    # f4(a_p,m,l,lp,mp,n,s)
    if a.y_axis:
        q1 = np.mgrid[0:rang_lp, 0:rang_s, 0:rang_n, 0:np.pi:rang_theta * 1j].astype('float32')
    else:
        q1 = np.mgrid[0:rang_lp, 0:rang_s, 0:rang_n, 0:rang_theta].astype('float32')
        q1[3, :] = q1[3, :] * np.pi / rang_theta

    q1[3, :] = np.cos(q1[3, :])

    q1[2, :] += a.n_min
    q1[2, :] = np.concatenate([q1[2, :, :, -a.n_min:, :],
                               q1[2, :, :, :-a.n_min, :]], axis=2)
    q1[1, :] += a.s_min
    q1[1, :] = np.concatenate([q1[1, :, -a.s_min:, :, :],
                               q1[1, :, :-a.s_min, :, :]], axis=1)

    q1a = bs.q_lsm(q1[0, :], q1[1, :], q1[2, :], q1[3, :])
    # del q1
    # q1(lp,s,n,theta)

    if a.y_axis:
        q2 = np.mgrid[0:rang_l, 0:rang_s, 0:rang_m, 0:np.pi:rang_theta * 1j].astype('float32')
    else:
        q2 = np.mgrid[0:rang_l, 0:rang_s, 0:rang_m, 0:rang_theta].astype('float32')
        q2[3, :] = q2[3, :] * np.pi / f3m.shape[2]

    q2[3, :] = np.cos(q2[3, :])

    q2[2, :] += a.m_min
    q2[2, :] = np.concatenate([q2[2, :, :, -a.m_min:, :],
                               q2[2, :, :, :-a.m_min, :]], axis=2)

    q2[1, :] += a.s_min
    q2[1, :] = np.concatenate([q2[1, :, -a.s_min:, :, :],
                               q2[1, :, :-a.s_min, :, :]], axis=1)

    q2a = bs.q_lsm(q2[0, :], q2[1, :], q2[2, :], q2[3, :])
    # del q2
    # q2(l,s,m,theta)

    '''measure or weight'''
    sintheta = np.sin(np.mgrid[0:np.pi:rang_theta * 1j])[np.newaxis, np.newaxis, :,
               np.newaxis, np.newaxis, np.newaxis,
               np.newaxis, np.newaxis]
    # sin_theta = np.sin(((np.mgrid[0:f3m.shape[2]]).astype('float') + 0.5) * np.pi / f3m.shape[2])

    '''sum to answer'''
    # standard(a_p, m, theta, lp, mp, n, s, l)
    # f3m(a_p, m, theta, lp, mp, n)
    f3b = f3m[:, :, :, :, :, :, np.newaxis, np.newaxis]
    # q1(lp,s,n,theta)
    q1b = q1a.transpose([3, 0, 2, 1])[np.newaxis, np.newaxis, :, :, np.newaxis, :, :, np.newaxis]

    # q2(l,s,m,theta)
    q2b = q2a.transpose([2, 3, 1, 0])[np.newaxis, :, :, np.newaxis, np.newaxis, np.newaxis, :, :]
    f_hat = np.zeros([rang_ap, rang_m, rang_lp, rang_mp, rang_s, rang_l]).astype('complex64')
    # for l in range(rang_l):
    #     for s in range(rang_s):
    #         for lp in range(rang_lp):
    #             f_hat[:, :, lp, :, s, l] = np.sum(f3b[:, :, :, lp, :, :, 0, 0] * q1b[:, :, :, lp, :, :, s, 0] *
    #                                               q2b[:, :, :, 0, :, :, s, l] * sintheta[:, :, :, 0, :, :, 0, 0],
    #                                               axis=(2, 4))
    for l in range(rang_l):
        for s in range(rang_s):
            for lp in range(rang_lp):
                f_hat[:, :, lp, :, s, l] = np.sum(np.trapz(f3b[:, :, :, lp, :, :, 0, 0] * q1b[:, :, :, lp, :, :, s, 0] *
                                                           q2b[:, :, :, 0, :, :, s, l] * sintheta[:, :, :, 0, :, :, 0,
                                                                                         0],
                                                           dx=np.pi / (rang_theta - 1), axis=2), axis=3)
    # f_hat *= np.pi / rang_theta
    # f4(a_p,m,lp,mp,n,s,l)
    # f_hat = np.sum(f4, axis=4)
    # f_hat(a_p,m,lp,mp,s,l)
    print 'Goodbye usphint. \n'
    return f_hat


def mltp(phi1, phi2):
    print 'This is mltp... \n'
    # phi(a_p,m,lp,mp,s,l)
    print phi1.shape
    # pdb.set_trace()
    m_mltp_max = min(a.m_max, a.mp_max)
    m_mltp_min = max(a.m_min, a.mp_min)
    rang_l_mltp = min(phi1.shape[2],phi1.shape[5])
    # rang_m = min(phi1.shape[1], phi1.shape[3])
    phi1 = phi1[:, :, :rang_l_mltp, :, :, :]
    phi1 = np.concatenate([phi1[:, :, :, :m_mltp_max + 1, :, :], phi1[:, :, :, m_mltp_min:, :, :]], axis=3)
    phi2 = phi2[:, :, :, :, :, :rang_l_mltp]
    phi2 = np.concatenate([phi2[:, :m_mltp_max + 1, :, :, :, :], phi2[:, m_mltp_min:, :, :, :, :]], axis=1)
    # phi2 = phi2[:, :rang_m, :rang_l, :rang_m, :, :rang_l]
    # pdb.set_trace()
    res = np.zeros(phi1.shape,dtype='complex64')
    for a_p in range(phi1.shape[0]):
        for s in range(phi1.shape[4]):
            res[a_p, :, :, :, s, :] = np.tensordot(phi2[a_p, :, :, :, s, :], phi1[a_p, :, :, :, s, :],
                                                   axes=([3, 0], [1, 2])).transpose([2, 0, 1, 3])
    print 'Goodbye mltp. \n'
    return res

@jit
def usphinti(phi):
    # phi(a_p,m,lp,mp,s,l)
    #     input shape expectations
    # a_p
    # m odd
    # lp (mp+1)/2
    # mp odd
    # s mp same
    # l (m+1)/2
    print 'This is ushpinti... \n'
    rang_ap, rang_m, rang_lp, rang_mp, rang_s, rang_l = phi.shape
    assert rang_m % 2 == 1 and rang_mp % 2 == 1 and rang_mp == rang_s == 2 * rang_lp - 1 \
           and rang_m == 2*rang_l-1

    '''calculation of basis function'''

    rang_theta = a.rang_theta
    rang_n = rang_mp
    a.n_imax = a.mp_max
    a.n_imin = a.mp_min
    rang_phi = a.rang_phi

    if a.y_axis:
        q1 = np.mgrid[0:rang_l, 0:rang_s, 0:rang_m, 0:np.pi:rang_theta * 1j].astype('float32')
    else:
        q1 = np.mgrid[0:rang_l, 0:rang_s, 0:rang_m, 0:rang_theta].astype('float32')
        q1[3, :] *= np.pi / rang_theta
    q1[3, :] = np.cos(q1[3, :])
    # TODO whoever cares the quadrature please add a proper weight.

    q1[2, :] += a.m_min
    q1[2, :] = np.concatenate([q1[2, :, :, -a.m_min:, :],
                               q1[2, :, :, :-a.m_min, :]], axis=2)
    q1[1, :] += a.s_min
    q1[1, :] = np.concatenate([q1[1, :, -a.s_min:, :, :],
                               q1[1, :, :-a.s_min, :, :]], axis=1)

    q1a = bs.q_lsm(q1[0, :], q1[1, :], q1[2, :], q1[3, :])
    # del q1
    # q1(l,s,m,theta)

    '''fft to g11'''
    # standard(a_p,m,lp,mp,s,l,theta)
    # phi(a_p,m,lp,mp,s,l)
    f = phi[:, :, :, :, :, :, np.newaxis]
    print f.shape, phi.shape

    # q1(l,s,m,theta)
    q1b = q1a.transpose([2, 1, 0, 3])[np.newaxis, :, np.newaxis, np.newaxis, :, :, :]
    g_full_m = f * q1b
    # del q1b
    # g_trunc_m = np.concatenate([g_full_m[:, :a.rang_phi / 2, :, :, :, :, :],
    #                             np.zeros([rang_ap, 1, rang_lp, rang_mp, rang_s,
    #                                       rang_l, rang_theta]),
    #                             g_full_m[:, -a.rang_phi / 2 + 1:, :, :, :, :, :]],
    #                            axis=1)
    g = np.concatenate([g_full_m[:,:a.m_max+1,:,:,:,:],np.zeros([rang_ap,rang_phi-rang_m,rang_lp,rang_mp,rang_s,rang_l,rang_theta]),
                        g_full_m[:,a.m_min:,:,:,:,:]],axis=1)
    g11 = ft.fft(g, axis=1)
    # del g_trunc_m
    # del g_full_m
    # g11(a_p,phi,lp,mp,s,l,theta)
    '''sum for g12'''
    # g12 = ne.evaluate('sum(g11, axis=5)')
    g12 = np.sum(g11,axis=5)
    # del g11
    # g12(a_p,phi,lp,mp,s,theta)

    '''multiply for g1'''
    g12a = g12[:, :, :, :, :, :, np.newaxis]
    # del g12
    # g12a(a_p,phi,lp,mp,s,theta,n)
    if a.y_axis:
        q2 = np.mgrid[0:rang_lp, 0:rang_s, 0:rang_n, 0:np.pi:rang_theta * 1j].astype('float32')
    else:
        q2 = np.mgrid[0:rang_lp, 0:rang_s, 0:rang_n, 0:rang_theta].astype('float32')
        q2[3, :] *= np.pi / rang_theta
        # q2(lp,s,n,theta)
    q2[3, :] = np.cos(q2[3, :])

    q2[2, :] += a.n_min
    q2[2, :] = np.concatenate([q2[2, :, :, -a.n_imin:, :],
                               q2[2, :, :, :-a.n_imin, :]], axis=2)
    q2[1, :] += a.s_min
    q2[1, :] = np.concatenate([q2[1, :, -a.s_min:, :, :],
                               q2[1, :, :-a.s_min, :, :]], axis=1)
    q2a = bs.q_lsm(q2[0, :], q2[1, :], q2[2, :], q2[3, :])
    # del q2
    # standard(a_p,phi,lp,mp,s,theta,n)
    # q2(lp,s,n,theta)
    q2b = q2a.transpose([0, 1, 3, 2])[np.newaxis, np.newaxis, :, np.newaxis, :, :, :]
    # del q2a

    ept = np.mgrid[0:rang_phi, 0:rang_n].astype('float32')
    # ept(phi,n)
    ept[0, :] *= 2 * np.pi / rang_phi
    ept[1, :] += a.n_imin
    ept[1, :] = np.concatenate([ept[1, :, -a.n_imin:],
                                ept[1, :, :-a.n_imin]], axis=1)
    epta = np.exp(1j * ept[0, :] * ept[1, :])

    eptb = epta[np.newaxis, :, np.newaxis, np.newaxis, np.newaxis, np.newaxis, :]
    # del epta

    g1 = g12a * q2b * eptb
    # g1(a_p,phi,lp,mp,s,theta,n)
    print 'Goodbye usphinti. \n'
    return g1.transpose([0, 5, 1, 2, 3, 4, 6])
    # (a_p,theta,phi,lp,mp,s,n)


# @profile
@jit
def sph2cart(phi):
    # g1(a_p,theta,phi,lp,mp,s,n)
    # input expectations:
    #    a_p
    #    theta odd
    #    phi doubly even
    #    lp, mp, s, n
    print 'This is sph2cart... \n'
    rang_ap, rang_theta, rang_phi, rang_lp, rang_mp, rang_s, rang_n = phi.shape
    f = np.concatenate([phi, phi[:, :, 0:1, :, :, :, :]], axis=2)
    # del phi
    assert rang_theta % 2 == 1 and rang_phi % 4 == 0
    rang_p = a.rang_p[0]
    print 'rang_p', rang_p
    print 'a.rang_p', a.rang_p
    p_max = (rang_p - 1) / 2
    p_min = -p_max
    axis_ticks = (np.mgrid[0:rang_ap], np.mgrid[0:np.pi:rang_theta * 1j], np.mgrid[0:2 * np.pi:(rang_phi + 1) * 1j])
    cart_grid = (np.mgrid[0:rang_p, 0:rang_p, 0:rang_p] + p_min).astype('float32')
    cart_grid += 0.
    p_x, p_y, p_z = cart_grid[0], cart_grid[1], cart_grid[2]
    x2y2 = p_x ** 2 + p_y ** 2
    #  print a.factor
    sph_grid = np.array([np.sqrt(x2y2 + p_z ** 2) / a.factor,
                         np.arctan2(np.sqrt(x2y2), p_z),
                         np.arctan2(p_y, p_x) % (2 * np.pi)]).transpose([1, 2, 3, 0])
    # print np.max(sph_grid[:,:,:,0])
    res = np.zeros([rang_p, rang_p, rang_p, rang_lp, rang_mp, rang_s, rang_n]).astype('complex64')
    for lp in range(rang_lp):
        for mp in range(rang_mp):
            for s in range(rang_s):
                for n in range(rang_n):
                    # print rang_p
                    # print rang_ap
                    res[:, :, :, lp, mp, s, n] = \
                        ip.interpn(axis_ticks, f[:, :, :, lp, mp, s, n], sph_grid)
    print 'Goodbye sph2cart'
    return res  # res(p1,p2,p3,lp,mp,s,n)


def e_ifft3(phi):
    print 'This is e_ifft3... \n'
    f = ft.ifftshift(phi, axes=(0, 1, 2))
    f = ft.fftn(f, axes=(0, 1, 2))
    f = np.concatenate([f[-a.org[0]:, :, :, :, :, :, :], f[:-a.org[0], :, :, :, :, :, :]], axis=0)
    f = np.concatenate([f[:, -a.org[1]:, :, :, :, :, :], f[:, :-a.org[1], :, :, :, :, :]], axis=1)
    f = np.concatenate([f[:, :, -a.org[2]:, :, :, :, :], f[:, :, :-a.org[2], :, :, :, :]], axis=2)
    print 'Goodbye e_ifft3. \n'
    return f


# @profile
# @profile
@jit
def so3inti(phi):  # phi(x1,x2,x3,lp,mp,s,n)
    #    expectation for inputs
    # x1,x2,x3 same, odd
    # lp (mp+1)/2
    # mp odd
    # s, mp, n same
    print 'This is so3inti... \n'
    f = phi
    rang_1, rang_2, rang_3, rang_lp, rang_mp, rang_s, rang_n = f.shape
    # p_lp,n,mp(cos_beta)
    # (x1,x2,x3,lp,mp,s,n,beta)
    rang_ga = rang_mp
    rang_al = rang_n
    rang_beta = a.rang_beta

    fa = f[:, :, :, :, :, :, :, np.newaxis]

    pa = np.mgrid[0:rang_lp, 0:rang_n, 0:rang_mp, 0:rang_beta].astype('float32')
    pa[3, :] += 0.5
    pa[3, :] *= np.pi / rang_beta
    pa[3, :] = np.cos(pa[3, :])
    pa[2, :] += a.mp_min
    pa[2, :] = np.concatenate([pa[2, :, :, -a.mp_min:, :],
                               pa[2, :, :, :-a.mp_min, :]], axis=2)
    pa[1, :] += a.n_imin
    pa[1, :] = np.concatenate([pa[1, :, -a.n_imin:, :, :],
                               pa[1, :, :-a.n_imin, :, :]], axis=1)

    # pb(lp, n, mp, beta)
    pb = bs.p_lmn(pa[0, :], pa[1, :], pa[2, :], pa[3, :])
    # del pa
    p = pb.transpose([0, 2, 1, 3])[np.newaxis, np.newaxis, np.newaxis, :, :, np.newaxis, :, :]
    res = np.zeros([rang_1, rang_2, rang_3, rang_ga, rang_al, rang_beta]).astype('complex64')
    for x1 in range(rang_1):
        for x2 in range(rang_2):
            for x3 in range(rang_3):
                f1 = p[0, 0, 0, :, :, :, :, :] * fa[x1, x2, x3, :, :, :, :, :]

                # fa *= p
                # f1 = fa
                f1[:, 1:a.mp_max + 1:2, :, :, :] *= -1
                f1[:, -1:a.mp_min - 1:-2, :, :, :] *= -1
                f1[:, :, :, 1:a.n_imax + 1:2, :] *= -1
                f1[:, :, :, -1:a.n_imin - 1:-2, :] *= -1
                # f1[:, :, :, :, 1:a.mp_max + 1:2, :, :, :] *= -1
                # f1[:, :, :, :, -1:a.mp_min - 1:-2, :, :, :] *= -1
                # f1[:, :, :, :, :, :, 1:a.n_imax + 1:2, :] *= -1
                # f1[:, :, :, :, :, :, -1:a.n_imin - 1:-2, :] *= -1

                # f2 = ft.fft(f1, axis=6)  # (x1,x2,x3,lp,mp,s,alpha,beta)
                f1 = ft.fft(f1, axis=3)  # (x1,x2,x3,lp,mp,s,alpha,beta)
                # del f1
    # fft of mp unchanged because structure of mp is unchanged.
                f1 = ft.fft(f1, axis=1)  # (x1,x2,x3,lp,gamma,s,alpha,beta)
                # f3 = ft.fft(f2, axis=4)  # (x1,x2,x3,lp,gamma,s,alpha,beta)
                # del f2
                f1 = np.sum(f1, axis=0)  # (x1,x2,x3,gamma,s,alpha,beta)
                # f4 = np.sum(f3, axis=3)  # (x1,x2,x3,gamma,s,alpha,beta)
                # del f3
                f1 = np.sum(f1, axis=1)  # (x1,x2,x3,gamma,alpha,beta)
                # f5 = np.sum(f4, axis=4)  # (x1,x2,x3,gamma,alpha,beta)
                # del f4
    # f5*=1/(2*np.pi**2)
                res[x1, x2, x3, :, :, :] = f1[:, :, :]

    f6 = res.transpose([0, 1, 2, 4, 5, 3])
    '''Real* Normalization'''
    f6 *= (2 * np.pi) ** 2
    print 'Goodbye so3inti. \n'
    return f6


def sph2cart3(phi):
    # g1(a_p,theta,phi,lp,mp,s,n)
    # input expectations:
    #    a_p
    #    theta odd
    #    phi doubly even
    #    lp, mp, s, n
    print 'This is sph2cart... \n'
    rang_ap, rang_theta, rang_phi, rang_lp, rang_mp, rang_s, rang_n = phi.shape
    f = np.concatenate([phi, phi[:, :, 0:1, :, :, :, :]], axis=2)
    f[:, 0, :, :, :, :, :] = f[:, 0, 0, np.newaxis, :, :, :, :]
    f[:, -1, :, :, :, :, :] = f[:, -1, 0, np.newaxis, :, :, :, :]
    f[0, :, :, :, :, :, :] = f[0, 0, 0, :, :, :, :]
    assert rang_theta % 2 == 1 and rang_phi % 4 == 0
    rang_p = a.rang_p[0]
    print 'rang_p', rang_p
    print 'a.rang_p', a.rang_p
    p_max = (rang_p - 1) / 2
    p_min = -p_max
    axis_ticks = (np.mgrid[0:rang_ap], np.mgrid[0:np.pi:rang_theta * 1j], np.mgrid[0:2 * np.pi:(rang_phi + 1) * 1j])
    cart_grid = (np.mgrid[0:rang_p, 0:rang_p, 0:rang_p] + p_min).astype('float64')
    cart_grid += 0.
    p_x, p_y, p_z = cart_grid[0], cart_grid[1], cart_grid[2]
    x2y2 = p_x ** 2 + p_y ** 2
    sph_grid = np.array([np.sqrt(x2y2 + p_z ** 2) / a.factor,
                         np.arctan2(np.sqrt(x2y2), p_z),
                         np.arctan2(p_y, p_x) % (2 * np.pi)]).transpose([1, 2, 3, 0])
    # print np.max(sph_grid[:,:,:,0])
    res = np.zeros([rang_p, rang_p, rang_p, rang_lp, rang_mp, rang_s, rang_n]).astype('complex64')
    for lp in range(rang_lp):
        for mp in range(rang_mp):
            for s in range(rang_s):
                for n in range(rang_n):
                    f_tmp = f[:, :, :, lp, mp, s, n]
                    f_tmp_r = np.real(f_tmp).astype('float64')
                    f_tmp_i = np.imag(f_tmp).astype('float64')
                    itp_r = tricubic.tricubic(list(f_tmp_r), [rang_ap, rang_theta, rang_phi + 1])
                    itp_i = tricubic.tricubic(list(f_tmp_i), [rang_ap, rang_theta, rang_phi + 1])
                    for p1 in range(rang_p):
                        for p2 in range(rang_p):
                            for p3 in range(rang_p):
                                res[p1, p2, p3, lp, mp, s, n] = itp_r.ip([sph_grid[p1, p2, p3, 0],
                                                                          sph_grid[p1, p2, p3, 1],
                                                                          sph_grid[p1, p2, p3, 2]]) + \
                                                                itp_i.ip([sph_grid[p1, p2, p3, 0],
                                                                          sph_grid[p1, p2, p3, 1],
                                                                          sph_grid[p1, p2, p3, 2]]) * 1j
                                # res[:, :, :, lp, mp, s, n] = \
                                #     ip.interpn(axis_ticks, f[:, :, :, lp, mp, s, n], sph_grid)
    print 'Goodbye sph2cart'
    return res  # res(p1,p2,p3,lp,mp,s,n)


# @profile
@jit
def inv_int(phi):  # integrated invert transformation wrapping in a loop of s,
    # reducing the memory consumption by 1 dimension.
    # phi(a_p,m,lp,mp,s,l)
    rang_ap, rang_m, rang_lp, rang_mp, rang_s, rang_l = phi.shape
    assert rang_m % 2 == 0 and rang_mp % 2 == 1 and rang_mp == rang_s == 2 * rang_lp - 1 \
           and rang_l == rang_lp
    rang_theta = a.rang_theta
    rang_n = rang_mp
    a.n_imax = a.mp_max
    a.n_imin = a.mp_min
    rang_phi = a.rang_phi
    rang_1, rang_2, rang_3 = a.rang_p[0], a.rang_p[1], a.rang_p[2]
    rang_ga = rang_mp
    rang_al = rang_n
    rang_beta = a.rang_beta
    res2 = np.zeros([rang_1, rang_2, rang_3, rang_ga, rang_al, rang_beta]).astype('complex64')
    for s in np.arange(rang_s) + a.s_min:

        '''    Previously usphinti    '''

        '''calculation of basis function'''

        if a.y_axis:
            q1 = np.mgrid[0:rang_l, 0:rang_m, 0:np.pi:rang_theta * 1j].astype('float32')
        q1[2, :] = np.cos(q1[2, :])
        # TODO whoever cares the quadrature please add a proper weight.
        q1[1, :] += a.m_min
        q1[1, :] = np.concatenate([q1[1, :, -a.m_min:, :],
                                   q1[1, :, :-a.m_min, :]], axis=1)
        qs1 = np.zeros(q1[1].shape).astype('float32')
        qs1[:] = s
        q1a = bs.q_lsm(q1[0, :], qs1, q1[1, :], q1[2, :])
        # del q1
        # q1(l,m,theta)

        '''fft to g11'''
        # standard(a_p,m,lp,mp,l,theta)
        # phi(a_p,m,lp,mp,s,l)
        f = phi[:, :, :, :, s, :, np.newaxis]
        print f.shape, phi.shape

        # q1(l,m,theta)
        q1b = q1a.transpose([1, 0, 2])[np.newaxis, :, np.newaxis, np.newaxis, :, :]
        g_full_m = f * q1b
        # del q1b
        g11 = ft.fft(g_full_m, axis=1)
        # g11(a_p,phi,lp,mp,l,theta)
        # del g_full_m
        '''sum for g12'''
        g12 = ne.evaluate('np.sum(g11, axis=4)')
        # del g11
        # g12(a_p,phi,lp,mp,theta)
        '''multiply for g1'''
        g12a = g12[:, :, :, :, :, np.newaxis]
        # g12a(a_p,phi,lp,mp,theta,n)
        if a.y_axis:
            q2 = np.mgrid[0:rang_lp, 0:rang_n, 0:np.pi:rang_theta * 1j].astype('float32')
            # q2(lp,n,theta)
        q2[2, :] = np.cos(q2[2, :])

        q2[1, :] += a.n_min
        q2[1, :] = np.concatenate([q2[1, :, -a.n_imin:, :],
                                   q2[1, :, :-a.n_imin, :]], axis=1)
        qs2 = np.zeros(q2[1].shape).astype('float32')
        qs2[:] = s
        q2a = bs.q_lsm(q2[0, :], qs2, q2[1, :], q2[2, :])
        # del q2
        # standard(a_p,phi,lp,mp,theta,n)
        # q2(lp,n,theta)
        q2b = q2a.transpose([0, 2, 1])[np.newaxis, np.newaxis, :, np.newaxis, :, :]
        # del q2a

        ept = np.mgrid[0:rang_phi, 0:rang_n].astype('float32')
        # ept(phi,n)
        ept[0, :] *= 2 * np.pi / rang_phi
        ept[1, :] += a.n_imin
        ept[1, :] = np.concatenate([ept[1, :, -a.n_imin:],
                                    ept[1, :, :-a.n_imin]], axis=1)
        epta = np.exp(1j * ept[0, :] * ept[1, :])

        eptb = epta[np.newaxis, :, np.newaxis, np.newaxis, np.newaxis, :]
        # del epta
        g1 = g12a * q2b * eptb
        # g1(a_p,phi,lp,mp,theta,n)
        g1 = g1.transpose([0, 4, 1, 2, 3, 5])
        '''    Previously sph2cart    '''
        # g1(a_p,theta,phi,lp,mp,n)
        f = np.concatenate([g1, g1[:, :, 0:1, :, :, :]], axis=2)
        rang_p = a.rang_p[0]
        print 'rang_p', rang_p
        print 'a.rang_p', a.rang_p
        p_max = (rang_p - 1) / 2
        p_min = -p_max
        axis_ticks = (np.mgrid[0:rang_ap], np.mgrid[0:np.pi:rang_theta * 1j], np.mgrid[0:2 * np.pi:(rang_phi + 1) * 1j])
        cart_grid = (np.mgrid[0:rang_p, 0:rang_p, 0:rang_p] + p_min).astype('float32')
        cart_grid += 0.
        p_x, p_y, p_z = cart_grid[0], cart_grid[1], cart_grid[2]
        x2y2 = p_x ** 2 + p_y ** 2
        print a.factor
        sph_grid = np.array([np.sqrt(x2y2 + p_z ** 2) / a.factor,
                             np.arctan2(np.sqrt(x2y2), p_z),
                             np.arctan2(p_y, p_x) % (2 * np.pi)]).transpose([1, 2, 3, 0])
        res = np.zeros([rang_p, rang_p, rang_p, rang_lp, rang_mp, rang_n]).astype('complex64')
        for lp in range(rang_lp):
            for mp in range(rang_mp):
                for n in range(rang_n):
                    # print rang_p
                    # print rang_ap
                    res[:, :, :, lp, mp, n] = \
                        ip.interpn(axis_ticks, f[:, :, :, lp, mp, n], sph_grid)
        # res(p1,p2,p3,lp,mp,n)
        # del g1, f
        '''    Previously e_ifft3    '''
        f = ft.ifftshift(res, axes=(0, 1, 2))
        f = ft.fftn(f, axes=(0, 1, 2))
        f = np.concatenate([f[-a.org[0]:, :, :, :, :, :], f[:-a.org[0], :, :, :, :, :]], axis=0)
        f = np.concatenate([f[:, -a.org[1]:, :, :, :, :], f[:, :-a.org[1], :, :, :, :]], axis=1)
        f = np.concatenate([f[:, :, -a.org[2]:, :, :, :], f[:, :, :-a.org[2], :, :, :]], axis=2)
        # f(x1,x2,x3,lp,mp,n)
        '''    Previously so3inti    '''

        # p_lp,n,mp(cos_beta)
        # (x1,x2,x3,lp,mp,n,beta)
        fa = f[:, :, :, :, :, :, np.newaxis]

        pa = np.mgrid[0:rang_lp, 0:rang_n, 0:rang_mp, 0:rang_beta].astype('float32')
        pa[3, :] += 0.5
        pa[3, :] *= np.pi / rang_beta
        pa[3, :] = np.cos(pa[3, :])
        pa[2, :] += a.mp_min
        pa[2, :] = np.concatenate([pa[2, :, :, -a.mp_min:, :],
                                   pa[2, :, :, :-a.mp_min, :]], axis=2)
        pa[1, :] += a.n_imin
        pa[1, :] = np.concatenate([pa[1, :, -a.n_imin:, :, :],
                                   pa[1, :, :-a.n_imin, :, :]], axis=1)
        # pb(lp, n, mp, beta)
        pb = bs.p_lmn(pa[0, :], pa[1, :], pa[2, :], pa[3, :])
        p = pb.transpose([0, 2, 1, 3])[np.newaxis, np.newaxis, np.newaxis, :, :, :, :]

        for x1 in range(rang_1):
            for x2 in range(rang_2):
                for x3 in range(rang_3):
                    # (lp,mp,n,beta)
                    f1 = p[0, 0, 0, :, :, :, :] * fa[x1, x2, x3, :, :, :, :]
                    f1[:, 1:a.mp_max + 1:2, :, :] *= -1
                    f1[:, -1:a.mp_min - 1:-2, :, :] *= -1
                    f1[:, :, 1:a.n_imax + 1:2, :] *= -1
                    f1[:, :, -1:a.n_imin - 1:-2, :] *= -1
                    f1 = ft.fft(f1, axis=2)  # (lp,mp,alpha,beta)
                    f1 = ft.fft(f1, axis=1)  # (lp,gamma,alpha,beta)
                    f1 = ne.evaluate('np.sum(f1, axis=0)')  # (gamma,alpha,beta)
                    res2[x1, x2, x3, :, :, :] += f1[:, :, :]
    res2 = res2.transpose([0, 1, 2, 4, 5, 3])
    '''Real * Normalization'''
    res2 *= (2 * np.pi) ** 2
    return res2
