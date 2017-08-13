#coding=UTF-8
"""Subroutine holding the calculation of basis functions"""

import numpy as np
import scipy.special as ss
import pdb
import warnings

warnings.filterwarnings("error", module='basis')


def p_lmn(l,mm,nn,z):
    m=np.copy(mm)
    n=np.copy(nn)
    pm = m+n < 0
    m[pm], n[pm] = -m[pm], -n[pm]
    cmp_m = m < n
    m[cmp_m], n[cmp_m] = n[cmp_m], m[cmp_m]
    if not (np.all(np.isfinite(m)) and np.all(np.isfinite(n))):
        print 'NaN before basis function is calculated!'
    tmp = np.zeros(l.shape)
    cmp_l = np.logical_and(l >= abs(m), l >= abs(n))
    tmp[cmp_l] = 2 ** (-m[cmp_l]) * np.sqrt(ss.factorial(l[cmp_l] - m[cmp_l]) * ss.factorial(l[cmp_l] + m[cmp_l]) /
                                            ss.factorial(l[cmp_l] - n[cmp_l]) / ss.factorial(l[cmp_l] + n[cmp_l])) * \
                 (1 - z[cmp_l]) ** ((m[cmp_l] - n[cmp_l]) / 2) * (1 + z[cmp_l]) ** ((m[cmp_l] + n[cmp_l]) / 2) * \
                 ss.eval_jacobi(l[cmp_l] - m[cmp_l], m[cmp_l] - n[cmp_l], m[cmp_l] + n[cmp_l], z[cmp_l])
    tmp[pm] *= (-1) ** (m[pm] - n[pm])
    tmp[cmp_m] *= (-1) ** (m[cmp_m] + n[cmp_m])
    # tmp[cmp_l] = np.where(np.logical_and(l>=abs(m),l>=abs(n)), np.where(pm, (-1)**(m-n), 1)*np.where(cmp,(-1)**(m+n),1)*2**(-m)*np.sqrt(ss.factorial(l-m)*ss.factorial(l+m)/ss.factorial(l-n)/ss.factorial(l+n))*(1-z)**((m-n)/2)*(1+z)**((m+n)/2)*ss.eval_jacobi(l-m,m-n,m+n,z), 0)

    if not np.all(np.isfinite(tmp)):
        print np.logical_not(np.isfinite(tmp))
        print 'NaN in normal indices!'
    return np.nan_to_num(tmp)


def q_lsm(l,s,m,cz):
    return (-1)**(l+s)*np.sqrt((2*l+1)/4/np.pi)*p_lmn(l,-s,m,cz)
    # return p_lmn(l,-s,m,cz)
