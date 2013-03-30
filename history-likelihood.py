"""
Estimate a likelihood from sampled histories.

This script should not compute the expm directly.
Input schema.
table wait
initial integer, final integer, state integer, wait real
table usage
initial integer, final integer, source integer, sink integer, usage real
"""

import argparse
import sqlite3
import math
import itertools

import numpy as np
import scipy.stats


#XXX this could go into a module of recipes
def pairwise(iterable):
    """
    This is an itertools recipe.
    s -> (s0,s1), (s1,s2), (s2, s3), ...
    """
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)

#XXX this should go into a separate module
def assert_stochastic_vector(v):
    if np.any(v < 0) or np.any(1 < v):
        raise Exception(
                'entries of a finite distribution vector should be in '
                'the inclusive interval [0, 1]')
    if not np.allclose(np.sum(v), 1):
        raise Exception(
                'entries of a finite distribution vector should sum to 1')

#XXX this should go into a separate module
def assert_rate_matrix(Q):
    if not np.allclose(np.sum(Q, axis=1), 0):
        raise Exception('expected rate matrix rows to sum to zero')
    if np.any(np.diag(Q) > 0):
        raise Exception('expected rate matrix diagonals to be non-positive')
    if np.any(Q - np.diag(np.diag(Q)) < 0):
        raise Exception('expected rate matrix off-diagonals to be non-negative')

#XXX this should go into a separate module
def assert_equilibrium(Q, distn):
    if not np.allclose(np.dot(distn, Q), 0):
        raise Exception('the distribution is not at equilibrium')

#XXX this should go into a separate module
def assert_detailed_balance(Q, distn):
    S = (Q.T * distn).T
    if not np.allclose(S, S.T):
        raise Exception('the detailed balance equations are not met')

#XXX this is copypasted
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


#XXX this is copypasted
def get_rate_matrix_info(cursor):
    """
    @param cursor: sqlite3 database cursor
    @return: sorted state list, stationary distribution, dense rate matrix
    """

    # get the sorted list of all states
    cursor.execute(
            'select state from distn '
            'union '
            'select source from rates '
            'union '
            'select sink from rates '
            'union '
            'select state from states '
            )
    states = sorted(t[0] for t in cursor)

    # count the states
    nstates = len(states)

    # define the map from state to rate matrix index
    s_to_i = dict((s, i) for i, s in enumerate(states))

    # construct the rate matrix
    cursor.execute('select source, sink, rate from rates')
    pre_Q = np.zeros((nstates, nstates), dtype=float)
    for si, sj, rate in cursor:
        pre_Q[s_to_i[si], s_to_i[sj]] = rate
    Q = pre_Q - np.diag(np.sum(pre_Q, axis=1))

    # construct the distribution of states at the root
    cursor.execute('select state, prob from distn')
    state_prob_pairs = list(cursor)
    distn = np.zeros(nstates)
    for state, prob in state_prob_pairs:
        distn[s_to_i[state]] = prob

    # assert that the distribution has the right form
    assert_stochastic_vector(distn)

    # assert that the rate matrix is actually a rate matrix
    assert_rate_matrix(Q)

    # assert that the distribution is at equilibrium w.r.t. the rate matrix
    assert_equilibrium(Q, distn)

    # assert that the detailed balance equations are met
    assert_detailed_balance(Q, distn)

    # return the validated inputs describing the stochastic process
    return states, distn, Q


def main(args):

    # read and validate the rate matrix info from the sqlite3 database file
    conn = sqlite3.connect(args.rates)
    cursor = conn.cursor()
    states, distn, Q = get_rate_matrix_info(cursor)
    conn.close()

    # read the path histories from the sqlite3 database file
    conn = sqlite3.connect(args.histories)
    cursor = conn.cursor()
    cursor.execute('select history, segment, state, blen from histories')
    data = sorted(cursor)
    cursor.execute('select state from histories')
    states = sorted(set(t[0] for t in cursor))
    conn.close()

    # summarize the Q matrix in a way that does not use expm
    rates, P = decompose_rates(Q)

    # Define the map from state to rate matrix index,
    # and count the number of different states in the history.
    s_to_i = dict((s, i) for i, s in enumerate(states))
    nstates = len(states)

    # parse the sample path histories from the table
    histories = []
    segment = []
    for row in data:
        history_index, segment_index, state, blen = row
        if history_index > len(histories):
            histories.append(segment)
            segment = []
        if history_index != len(histories):
            print row
            raise Exception('invalid history index')
        if segment_index != len(segment):
            raise Exception('invalid segment index')
        segment.append((state, blen))
    histories.append(segment)

    # compute the likelihood for each path history
    likelihoods = []
    for history in histories:
        likelihood = get_conditional_likelihood_from_history(
                Q, rates, P, states, history)
        likelihoods.append(likelihood)

    # report likelihood summary
    print 'likelihood arithmetic mean:', np.mean(likelihoods)
    print 'likelihood geometric mean:', scipy.stats.gmean(likelihoods)
    print 'likelihood harmonic mean:', scipy.stats.hmean(likelihoods)
    print


def get_conditional_likelihood_from_history(Q, rates, P, states, history):
    """
    Each history is a sample path.
    Each sample path is a sequence of (state, wait) pairs.
    Return the likelihood conditional on the initial state.
    @param Q: rate matrix
    @param rates: rate away from each state
    @param P: instantaneous transition probability distributions
    @param states: ordered states
    @param history: a single sample path
    @return: likelihood conditional on initial state
    """

    # FIXME: eventually do log likelihood instead

    # init the log likelihood
    likelihood = 1.0

    # Precompute the map from states to state indices.
    s_to_i = dict((s, i) for i, s in enumerate(states))
    n = len(states)

    # Add dwell time contributions.
    for state, wait in history:
        i = s_to_i[state]
        q = rates[i]
        likelihood *= math.exp(-q*wait)

    # Add transition contributions.
    for ((state_a, wait_a), (state_b, wait_b)) in pairwise(history):
        i = s_to_i[state_a]
        j = s_to_i[state_b]
        q = rates[i]
        likelihood *= q * P[i, j]

    # Return the likelihood.
    return likelihood


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rates', default='rate.matrix.db',
            help='time-reversible rate matrix as an sqlite3 database file')
    parser.add_argument('--histories', default='histories.db',
            help='input path histories as an sqlite3 database file')
    parser.add_argument('--outfile', default='averages.db',
            help='output averages as an sqlite3 database file')
    main(parser.parse_args())

