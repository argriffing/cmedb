"""
Compute endpoint-conditioned expectations.

The meanings of the components of the script name are as follows.
This is related to continuous-time Markov chain processes hence ctmc.
The processes are assumed to be statistically time-reversible hence rev.
The ec means endpoint-conditioned,
which is sometimes also referred to as bridged, like a Brownian bridge
or bridge sampling.
The process is assumed to act in continuous time along a path rather
than along a star tree or along a general tree,
hence the word 'path' in the script name.
Finally, we want to compute conditional expectations of some stuff.
The conditional expectations are the quantities in equation (4)
of Holmes and Rubin (2002).
These are the expected number of paths that start in state i,
the expected wait in state i (i.e. the amount of time spent in i)
and the expected usage of transition i->j.
"""

def get_decay_mode_interactions(w, T):
    """
    @param w: eigenvalues
    @param T: elapsed time
    @return: a matrix of interactions between pairs of decay modes
    """
    nstates = w.shape[0]
    J = np.empty((nstates, nstates))
    v = np.exp(w * T)
    for k in range(nstates):
        for l in range(nstates):
            denom = w[k] - w[l]
            if denom:
                J[k, l] = (v[k] - v[l]) / denom
            else:
                J[k, l] = T * v[k]
    return J

def foo(a, b, V, J):
    """
    I am not sure exactly what this is,
    @param a: initial endpoint state
    @param b: final endpoint state
    @param V: orthogonal matrix as an ndarray
    @param J: decay mode interactions including eigenvalue information
    @return: an expectation matrix
    """
    nstates = V.shape[0]
    T = np.empty_like(J)
    #XXX
    pass


