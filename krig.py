# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 19:41:26 2015

@author: yuhao
"""
from __future__ import division
import numpy as np
from pandas import DataFrame
from scipy.spatial.distance import pdist, squareform
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import matplotlib
import time
# from functools import partial

# import multiprocessing as mp
# from numba import jit


def speed_it(func):
    def wrapped(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        print "time elaped: {}".format(time.time() - start)
        return result
    return wrapped


def SVh(P, h, bw):
    '''
    Experimental semivariogram for a single lag
    '''
    pd = squareform(pdist(P[:, :2]))
    N = pd.shape[0]
    Z = list()
    for i in xrange(N):
        for j in xrange(i+1, N):
            if pd[i, j] >= h - bw and pd[i, j] <= h + bw:
                Z.append((P[i, 2] - P[j, 2])**2.0)
    return np.sum(Z) / (2.0 * len(Z))


def SV(P, hs, bw):
    '''
    Experimental variogram for a collection of lags
    '''
    sv = list()
    for h in hs:
        sv.append(SVh(P, h, bw))
    sv = [[hs[i], sv[i]] for i in range(len(hs)) if sv[i] > 0]
    return np.array(sv).T


def C(P, h, bw):
    """
    calculate nugget covariance of the variogram

    Parameters
    ----------
    P : ndarray
      samples
    h : scalar
      lag
    bw : scalar
      bandwidth

    Returns
    -------
    c0 : number
      nugget variance
    """
    c0 = np.var(P[:, 2])  # sill variance, which is a priori variance
    if h == 0:
        return c0
    return c0 - SVh(P, h, bw)


@speed_it
def opt(fct, x, y, c0, parameterRange=None, meshSize=1000):
    """
    Parameters
    ----------
    fct : callable
      functions to minimize
    x : ndarray
      input x array
    y : ndarray
      input y array
    c0 : scalar
      covariance

    Returns
    -------
    a : number
      the best fitting coefficient
    """
    if parameterRange is None:
        parameterRange = [x[1], x[-1]]

    def residuals(pp, yy, xx):
        a = pp
        err = yy - fct(xx, a, c0)
        return err

    p0 = x[1]
    plsq = leastsq(residuals, p0, args=(y, x))
    return plsq[0][0]


def spherical(h, a, c0):
    if h <= a:
        return c0 * (1.5*h/a - 0.5*(h/a)**3.0)
    else:
        return c0


def cvmodel(P, model, hs, bw):
    '''
    Parameters
    ----------
    P : ndarray
      data
    model : callable
        modeling function

        - spherical
        - exponential
        - gaussian
    hs : ndarray
      distances
    bw : scalar
      bandwidth

    Returns
    -------
    covfct : callable
      the optimized function modeling the semivariogram
    '''
    sv = SV(P, hs, bw)  # calculate the semivariogram
    C0 = C(P, hs[0], bw)  # calculate the nugget
    param = opt(model, sv[0], sv[1], C0)

    def covfct(h, a=param):
        return model(h, a, C0)
#    covfct = partial(model, a=param, c0=C0)
    return covfct


def krige(P, model, hs, bw, u, N):
    '''
    Parameters
    ----------
    P : ndarray
        data
    model : callable
        modeling function
        (spherical, exponential, gaussian)
    hs:
        kriging distances
    bw:
        kriging bandwidth
    u:
        unsampled point
    N:
        number of neighboring points to consider

    Returns
    -------
    estimate : scalar
        krigged value
    '''
#    covfct = cvmodel(P, model, hs, bw)  # covariance function
    covfct = model
    mu = np.mean(P[:, 2])  # mean of the variable
    # distance between u and each data point in P
    d = np.sqrt((P[:, 0] - u[0])**2.0 + (P[:, 1] - u[1])**2.0)
    P = np.vstack((P.T, d)).T  # add these distances to P
    # sort P by these distances
    # take the first N of them
    P = P[d.argsort()[:N]]
    k = covfct(P[:, 3])  # apply the covariance model to the distances
    k = np.matrix(k).T  # cast as a matrix
    # form a matrix of distances between existing data points
    K = squareform(pdist(P[:, :2]))
    # apply the covariance model to these distances
    K = covfct(K.ravel())
    K = np.array(K)
    K = K.reshape(N, N)
    K = np.matrix(K)
    # calculate the kriging weights
    weights = np.linalg.inv(K) * k
    weights = np.array(weights)

    variance_sample = np.var(P[:, 2])
    variance = variance_sample - np.dot(k.T, weights)

    variance = np.abs(variance)

    residuals = P[:, 2] - mu  # calculate the residuals

    # calculate the estimation
    estimation = np.dot(weights.T, residuals) + mu

    # vc = np.sqrt(variance) / estimation

    return float(estimation), float(variance)  # float(vc)

if __name__ == '__main__':
    # import and munge data
    z = open('WGTutorial/ZoneA.dat', 'r').readlines()
    z = [i.strip().split() for i in z[10:]]
    z = np.array(z, dtype=np.float)
    z = DataFrame(z, columns=['x', 'y', 'thk', 'por', 'perm',
                              'lperm', 'lpermp', 'lpermr'])

    P = np.array(z[['x', 'y', 'por']])
    bw = 500  # bandwidth, plus or minus 250 meters
    # lags in 500 meter increments from zero to 10,000
    hs = np.arange(0, 10500, bw)
    sv = SV(P, hs, bw)

#    fig1, ax1 = plt.subplots()
#    ax1.plot(sv[0], sv[1], '.-')
#    ax1.xlabel('Lag [m]')
#    ax1.ylabel('Semivariance')
#    ax1.title('Sample Semivariogram')
#    fig1.savefig('sample_semivariogram.png', fmt='png', dpi=200)

    vfunc_model = np.vectorize(spherical)
    sp = cvmodel(P, model=vfunc_model, hs=np.arange(0, 10500, 500), bw=500)

#    fig2, ax2 = plt.subplots()
#    ax2.plot(sv[0], sv[1], '.-')
#    ax2.plot(sv[0], sp(sv[0]))
#    ax2.set_title('Spherical Model')
#    ax2.set_ylabel('Semivariance')
#    ax2.set_xlabel('Lag [m]')
#    fig2.savefig('semivariogram_model.png', fmt='png', dpi=200)

    X0, X1 = P[:, 0].min(), P[:, 0].max()
    Y0, Y1 = P[:, 1].min(), P[:, 1].max()
    Z = np.zeros((80, 100))
    V = np.zeros((80, 100))
    dx, dy = (X1 - X0)/100.0, (Y1 - Y0)/80.0

    model = cvmodel(P, vfunc_model, hs, bw)

    start = time.time()
    m = 80
    for i in xrange(m):
        percent = int(i / m * 100)
        print '[' + percent * '=' + (100-percent) * ' ' + '] ' +\
            str(percent) + '%'
        for j in xrange(100):
            Z[i, j], V[i, j] = krige(P, model, hs, bw, (dy*j, dx*i), 16)
    print '[' + 100 * '=' + '] ' + '100%' + '\nCompleted!'
    #        Z[i, j]=answer[0]
    #        V[i, j]=answer[1]
    print "\nTotal time: ", (time.time() - start)

    # Plot
    matplotlib.use('Agg')
    cdict = {'red':   ((0.0, 1.0, 1.0),
                       (0.5, 225/255., 225/255.),
                       (0.75, 0.141, 0.141),
                       (1.0, 0.0, 0.0)),
             'green': ((0.0, 1.0, 1.0),
                       (0.5, 57/255., 57/255.),
                       (0.75, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),
             'blue':  ((0.0, 0.376, 0.376),
                       (0.5, 198/255., 198/255.),
                       (0.75, 1.0, 1.0),
                       (1.0, 0.0, 0.0))}
    my_cmap = matplotlib.colors.LinearSegmentedColormap(
                                                'my_colormap', cdict, 256)
#    fig, ax = plt.subplots()
#    H = np.zeros_like(Z)
#    for i in range(Z.shape[0]):
#        for j in range(Z.shape[1]):
#            H[i, j] = np.round(Z[i, j] * 3)
#    ax.matshow(H, cmap=my_cmap, interpolation='nearest')
#    ax.scatter(z.x/200.0, z.y/200.0, facecolor='none', linewidths=0.75, s=50)
#    ax.set_xlim(0, 99)
#    ax.set_ylim(0, 80)
#    ax.set_xticks([25, 50, 75], [5000, 10000, 15000])
#    ax.set_yticks([25, 50, 75], [5000, 10000, 15000])
#    fig.savefig('krigingpurple.png', fmt='png', dpi=200)
#
#    figx, axx = plt.subplots()
#    HV = np.zeros_like(V)
#    for i in range(V.shape[0]):
#        for j in range(V.shape[1]):
#            HV[i, j] = np.round(V[i, j] * 3)
#    axx.matshow(HV, cmap=my_cmap, interpolation='nearest')
#    axx.scatter(z.x/200.0, z.y/200.0, facecolor='none', linewidths=0.75, s=50)
#    axx.set_xlim(0, 99)
#    axx.set_ylim(0, 80)
#    axx.set_xticks([25, 50, 75], [5000, 10000, 15000])
#    axx.set_yticks([25, 50, 75], [5000, 10000, 15000])
#    figx.savefig('krigingpurple2.png', fmt='png', dpi=200)

    fig, ax = plt.subplots()
    H = np.zeros_like(Z)
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            H[i, j] = np.round(Z[i, j] * 3)
    ax.matshow(H, cmap=my_cmap, interpolation='nearest')
    ax.scatter(z.x/200.0, z.y/200.0, facecolor='none', linewidths=0.75, s=50)
    ax.set_xlim(0, 99)
    ax.set_ylim(0, 80)
    ax.set_xticks([25, 50, 75], [5000, 10000, 15000])
    ax.set_yticks([25, 50, 75], [5000, 10000, 15000])
    fig.savefig('krigingpurple.svg', fmt='svg', dpi=200)
