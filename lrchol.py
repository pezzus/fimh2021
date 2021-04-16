import numpy as np

def lowRankCholesky(dA,rowfun,M=200,tol=1e-4,verbose=False):
    """Compute low-rank pivoted Cholesky factorization.
    :param dA:     initial diagonal
    :param rowfun: function returning row i of matrix A
    :param M:      initial guess of reduced rank (M << N)
    :param tol:    tolerance
    """

    # diagonal
    d = dA.copy()

    # dimension of the matrix
    N = d.size

    # Guess reduced rank
    M = min(N,M)

    # Allocate space for efficiency
    L = np.zeros((M,N))

    # permutation list for pivoting
    p = np.arange(0,N)

    # initial error
    error = np.linalg.norm(d,1)

    # index from 0 and not from 1
    m = 0

    while error > tol:
        i = m + np.argmax( d[p[m:N]] )

        p[[m,i]] = p[[i,m]]

        L[m,p[m]] = np.sqrt(d[p[m]])

        # evaluate row p[m] of the full matrix
        a = rowfun(p[m])
        s = np.dot(L[0:m,p[m]],L[0:m,p[m+1:N]])
        L[m,p[m+1:N]] = (a[p[m+1:N]] - s)/L[m,p[m]]
        d[p[m+1:N]] -= L[m,p[m+1:N]]**2

        error = d[p[m+1:N]].sum()
        if verbose: print(m,error)

        if m+1 == M: break
        m += 1

    # we transpose so that the matrix is N x M
    L = np.delete(L,range(m,M),axis=0).T

    return L,m

