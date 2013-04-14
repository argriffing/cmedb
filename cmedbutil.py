"""
This is a temporary module for miscellaneous reusable functions.
"""

import itertools
import argparse

import numpy as np
import networkx as nx


###############################################################################
# Standard itertools recipes.

def pairwise(iterable):
    """
    This is an itertools recipe.
    s -> (s0,s1), (s1,s2), (s2, s3), ...
    """
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)


###############################################################################
# Not itertools recipes.


def first_element(it):
    for x in it:
        return x

def nth_element(n, seq):
    return seq[n]


###############################################################################
# Extra types for argparse.

def nonneg_int(x):
    x = int(x)
    if x < 0:
        raise argparse.ArgumentTypeError(
                'value must be a non-negative integer')
    return x

def pos_int(x):
    x = int(x)
    if x < 1:
        raise argparse.ArgumentTypeError(
                'value must be a positive integer')
    return x

def nonneg_float(x):
    x = float(x)
    if x < 0:
        raise argparse.ArgumentTypeError(
                'value must be a non-negative floating point number')
    return x

def pos_float(x):
    x = float(x)
    if x <= 0:
        raise argparse.ArgumentTypeError(
                'value must be a positive floating point number')
    return x

def prob_float(x):
    x = float(x)
    if not (0 <= x <= 1):
        raise argparse.ArgumentTypeError(
                'value must be in the closed interval [0, 1]')
    return x


###############################################################################
# Validation of invariants related to stochastic processes.

def assert_stochastic_vector(v):
    if np.any(v < 0) or np.any(1 < v):
        raise Exception(
                'entries of a finite distribution vector should be in '
                'the inclusive interval [0, 1] '
                'min: %s  max: %s' % (min(v), max(v)))
    if not np.allclose(np.sum(v), 1):
        raise Exception(
                'entries of a finite distribution vector should sum to 1')

def assert_rate_matrix(Q):
    if not np.allclose(np.sum(Q, axis=1), 0):
        raise Exception('expected rate matrix rows to sum to zero')
    if np.any(np.diag(Q) > 0):
        raise Exception('expected rate matrix diagonals to be non-positive')
    if np.any(Q - np.diag(np.diag(Q)) < 0):
        raise Exception('expected rate matrix off-diagonals to be non-negative')

def assert_equilibrium(Q, distn):
    if not np.allclose(np.dot(distn, Q), 0):
        raise Exception('the distribution is not at equilibrium')

def assert_detailed_balance(Q, distn):
    S = (Q.T * distn).T
    if not np.allclose(S, S.T):
        raise Exception('the detailed balance equations are not met')


###############################################################################
# Validation of invariants related to unrooted trees with edge lengths.

def assert_connected_acyclic_graph(G):
    """
    Check that the graph is connected and has no cycles.
    @param G: networkx undirected graph
    """
    if not nx.is_connected(G):
        raise Exception('the graph is not connected')
    if nx.cycle_basis(G):
        raise Exception('the graph has a cycle')


###############################################################################
# Extra functions of possibly general interest.

def random_category(distn):
    """
    Sample from a categorical distribution.
    Note that this is not the same as random.choice(distn).
    Maybe a function like this will eventually appear
    in python or numpy or scipy.
    @param distn: categorical distribution as a stochastic vector
    @return: category index as a python integer
    """
    nstates = len(distn)
    np_index = np.dot(np.arange(nstates), np.random.multinomial(1, distn))
    return int(np_index)


###############################################################################
# Functions specific to continuous time Markov processes.

def decompose_rates(Q):
    """
    Break a rate matrix into two parts.
    The first part consists of the rates away from each state;
    this information is contained in the diagonal of the rate matrix.
    The second part consists of a transition matrix
    that defines the distribution over sink states conditional
    on an instantaneous change away from a given source state.
    Note that this function never requires expm of Q.
    Also, P preserves the sparsity pattern of Q.
    @param Q: rate matrix
    @return: rates, P
    """
    nstates = len(Q)
    rates = -np.diag(Q)
    P = np.array(Q)
    for i, rate in enumerate(rates):
        if rate:
            P[i, i] = 0
            P[i] /= rate
    return rates, P



#FIXME: add testing

