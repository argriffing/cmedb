"""
Compute blinking process path history likelihood.

The process is assumed to consist of two simultaneously evolving components.
The primary component is assumed to be observable,
and the marginal process of this component
is a Markov modulated continuous-time process.
The secondary component is assumed to not be observable,
and it consists of multiple Markov-modulated marginal processes.
This script tries to compute the history likelihood
by integrating over the unobserved process.
This will involve using dynamic programming to compute
the likelihood for a discrete-time finite-state
inhomogeneous hidden Markov model.
.
Compute a partially observed path history likelihood.
This is a complicated script involving Markov modulation.
An observable genetic history is evolving along a path,
and we assume that we know its state at every point along this path.
In other words,
we know its state at the beginning of the path and
we know each time that its state changes.
If this observable genetic history were to have evolved according to a
continuous-time Markov process, then computing the path history likelihood
would be fairly straightforward.
On the other hand,
in this script we assume that the observable genetic history evolves
according to a Markov modulated continuous-time process
for which the modulating states are hidden.
.
The following command can be used to sample codon.path.history.db.
$ python ctmc-segment-bridge-sampling.py
--rates=codon.rate.matrix.db --method=modified-rejection
--outfile=codon.path.history.db --nsamples=1 --table=history --elapsed=10
--initial=0 --final=60
"""

import argparse
import sqlite3
import math
import itertools
from itertools import permutations

import numpy as np
import networkx as nx
import scipy.linalg

import cmedbutil


def get_micro_rate_matrix(rate_off, rate_on, rate_absorb):
    """
    The micro state order is ('off', 'on', 'absorbed').
    @param rate_off: rate from the 'on' state to the 'off' state
    @param rate_on: rate from the 'off' state to the 'on' state
    @param rate_absorb: rate from the 'on' state to the 'absorbing' state
    @return: a continuous time rate matrix
    """
    return np.array([
        [-rate_on, rate_on, 0],
        [rate_off, -(rate_off + rate_absorb), rate_absorb],
        [0, 0, 0]], dtype=float)


#XXX copypasted
def get_sparse_rate_matrix_info(cursor):
    """
    This is a non-numpy customization of the usual get_rate_matrix_info.
    In this function we ignore the 'states' table with the state names,
    but we care about the sparse rates and the equilibrium distribution.
    Return a couple of things.
    The first thing is a map from a state to an equilibrium probability.
    The second thing is a sparse rate matrix as a networkx weighted digraph.
    @param cursor: sqlite3 database cursor
    @return: eq distn, rate graph
    """

    # get a set of all states
    cursor.execute(
            'select state from distn '
            'union '
            'select source from rates '
            'union '
            'select sink from rates '
            'union '
            'select state from states '
            )
    states = set(t[0] for t in cursor)
    nstates = len(states)

    # get the sparse equilibrium distribution
    cursor.execute('select state, prob from distn')
    distn = dict(cursor)

    # construct the rate matrix as a networkx weighted directed graph
    dg = nx.DiGraph()
    cursor.execute('select source, sink, rate from rates')
    for a, b, weight in cursor:
        dg.add_edge(a, b, weight=weight)

    # assert that the distribution has the right form
    if not all(0 <= p <= 1 for p in distn.values()):
        raise Exception(
                'equilibrium probabilities '
                'should be in the interval [0, 1]')
    if not np.allclose(sum(distn.values()), 1):
        raise Exception('equilibrium probabilities should sum to 1')

    # assert that rates are not negative
    if any(data['weight'] < 0 for a, b, data in dg.edges(data=True)):
        raise Exception('rates should be non-negative')

    # assert detailed balance
    for a, b in permutations(states, 2):
        if b in dg[a] and a in dg[b]:
            if not np.allclose(
                    distn[a] * dg[a][b]['weight'],
                    distn[b] * dg[b][a]['weight'],
                    ):
                raise Exception('expected detailed balance')

    # return the eq distn and the rate graph
    return distn, dg


def get_blink_thread_ll(
        part, partition, distn, dg, path_history,
        rate_on, rate_off):
    """
    @param part: the part of the partition defining the current blink thred
    @param partition: a map from primary state to part
    @param distn: map from primary state to equilibrium probability
    @param dg: sparse primary state rate matrix as weighted directed networkx
    @param path_history: triples of (segment, state, blen)
    @param rate_on: a blink rate
    @param rate_off: a blink rate
    @return: a log likelihood
    """
    return -3.14


def main(args):

    # read the sparse rate matrix from a database file
    conn = sqlite3.connect(args.rates)
    cursor = conn.cursor()
    distn, dg = get_sparse_rate_matrix_info(cursor)
    conn.close()

    # Read the partition from a database file.
    # The hidden blinking process controls the state transition
    # of the primary process according to this partition
    # of primary process states.
    conn = sqlite3.connect(args.partition)
    cursor = conn.cursor()
    partition = dict(cursor.execute('select state, part from partition'))
    conn.close()

    # get the set of indices of parts
    parts = set(partition.values())
    nparts = len(parts)

    # read the primary state path history from the sqlite3 database file
    conn = sqlite3.connect(args.path_history)
    cursor = conn.cursor()
    cursor.execute('select segment, state, blen from history order by segment')
    path_history = list(cursor)
    conn.close()

    # compute the log likelihood by summing over blink thread log likelihoods
    ll = 0.0
    for part in range(nparts):
        ll += get_blink_thread_ll(
                part, partition, distn, dg, path_history,
                args.rate_on, args.rate_off)

    # report the likelihood
    print 'log likelihood:', ll
    print 'likelihood:', math.exp(ll)


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
    for ((state_a, wait_a), (state_b, wait_b)) in cmedbutil.pairwise(history):
        i = s_to_i[state_a]
        j = s_to_i[state_b]
        q = rates[i]
        likelihood *= q * P[i, j]

    # Return the likelihood.
    return likelihood


def get_conditional_final(Q, rates, P, states, history):
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

    # Precompute the map from states to state indices.
    s_to_i = dict((s, i) for i, s in enumerate(states))

    # Return the conditional likelihood.
    if len(history) == 1:
        return 1.0
    else:
        (state_a, wait_a), (state_b, wait_b) = history[-2:]
        i = s_to_i[state_a]
        j = s_to_i[state_b]
        return P[i, j]


if __name__ == '__main__':

    # define some methods
    method_choices = (
            'brute',
            'dynamic',
            )

    # define the command line parameters
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rate-on',
            type=cmedbutil.pos_float, default=1.0,
            help='rate at which blink states change from off to on')
    parser.add_argument('--rate-off',
            type=cmedbutil.pos_float, default=1.0,
            help='rate at which blink states change from on to off')
    parser.add_argument('--method', choices=method_choices, default='brute',
            help='method of integrating over hidden blink states')
    parser.add_argument('--path-history', default='path.history.db',
            help='input path history as an sqlite3 database file')
    parser.add_argument('--rates', default='rate.matrix.db',
            help=('continuous-time Markov substitution process '
                'as an sqlite3 database file'))
    parser.add_argument('--partition', default='partition.db',
            help=('a partition of the primary states '
                'as an sqlite3 database file'))
    main(parser.parse_args())

