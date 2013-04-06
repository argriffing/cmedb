"""
Compute slow blinking process path history likelihood.

This is a slow-blinking limit,
but the blinked-on proportion is still meaningfully defined.
"""

import argparse
import sqlite3
import math
import itertools
from itertools import permutations
from itertools import product

import numpy as np
import networkx as nx
import scipy.linalg

import cmedbutil


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


def get_partially_observed_blink_thread_log_likelihood(
        part, partition, distn, dg, path_history,
        proportion_on,
        blink_state,
        ):
    """
    @param part: the part of the partition defining the current blink thread
    @param partition: a map from primary state to part
    @param distn: map from primary state to equilibrium probability
    @param dg: sparse primary state rate matrix as weighted directed networkx
    @param path_history: triples of (segment, state, blen)
    @param proportion_on: equilibrium probability of blink state 'on'
    @param blink_state: the state of the blink thread
    @return: a log likelihood
    """

    # count the number of segments
    nsegments = len(path_history)

    # init the path log likelihood
    finite_path_ll = 0.0

    # slowly and inefficiently get the first primary state
    seg_seq, primary_state_seq, duration_seq = zip(*path_history)
    initial_primary_state = primary_state_seq[0]

    # The initial state contributes to the likelihood.
    if partition[initial_primary_state] == part:
        if not blink_state:
            return -float('Inf')
    else:
        blink_distn = np.array([1-proportion_on, proportion_on], dtype=float)
        log_blink_distn = np.log(blink_distn)
        finite_path_ll += log_blink_distn[blink_state]

    # multiply the likelihood across all segments along the thread
    for i in range(nsegments):

        # extract the primary and blink state
        segment, primary_state, duration = path_history[i]

        # If the blink state is on,
        # then penalize long duration in primary states with high escape rate.
        # Otherwise terminally penalize impossible primary states.
        if blink_state:
            rate_absorb = 0.0
            for sink in dg.successors(primary_state):
                rate = dg[primary_state][sink]['weight']
                if partition[sink] == part:
                    rate_absorb += rate
            finite_path_ll -= duration * rate_absorb
        else:
            if partition[primary_state] == part:
                return -float('Inf')

    # Return the finite path likelihood for this blinking thread.
    return finite_path_ll


def get_blink_thread_log_likelihood(
        part, partition, distn, dg, path_history,
        proportion_on,
        ):
    """
    @param part: the part of the partition defining the current blink thred
    @param partition: a map from primary state to part
    @param distn: map from primary state to equilibrium probability
    @param dg: sparse primary state rate matrix as weighted directed networkx
    @param path_history: triples of (segment, state, blen)
    @param proportion_on: equilibrium probability of blink state 'on'
    @return: log likelihood
    """

    # count the number of segments and the number of segment endpoints
    nsegments = len(path_history)
    npoints = nsegments + 1

    # initialize the likelihood
    log_likelihoods = []

    # sum the likelihood over all possible states at the segment endpoints
    for blink_state in (0, 1):
        log_likelihood = get_partially_observed_blink_thread_log_likelihood(
                part, partition, distn, dg, path_history,
                proportion_on,
                blink_state)
        log_likelihoods.append(log_likelihood)

    # report the log likelihood
    return scipy.misc.logsumexp(log_likelihoods)


def get_primary_log_likelihood(distn, dg, path_history):
    """
    This includes the initial state and the transitions.
    It does not include the dwell times,
    because these contributions are added by the blinking threads.
    @param distn: map from primary state to equilibrium probability
    @param dg: sparse primary state rate matrix as weighted directed networkx
    @param path_history: triples of (segment, state, blen)
    @return: log likelihood
    """

    # initialize
    log_likelihood = 0.0
    seg_seq, state_seq, duration_seq = zip(*path_history)

    # add the contribution of the equilibrium distribution
    initial_primary_state = state_seq[0]
    log_likelihood += math.log(distn[initial_primary_state])

    # add the contribution of the primary state transitions
    for a, b in cmedbutil.pairwise(state_seq):
        rate = dg[a][b]['weight']
        log_likelihood += math.log(rate)

    # return the log likelihood
    return log_likelihood


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

    # init the log likelihood
    log_likelihood = 0.0

    # add primary state log likelihood contributions
    log_likelihood += get_primary_log_likelihood(distn, dg, path_history)

    # add log likelihood blink thread log likelihood contributions
    for part in range(nparts):
        log_likelihood += get_blink_thread_log_likelihood(
                part, partition, distn, dg, path_history,
                args.proportion_on)

    # report the likelihood
    print 'log likelihood:', log_likelihood
    print 'likelihood:', math.exp(log_likelihood)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--proportion-on',
            type=cmedbutil.prob_float,
            help='on/(on+off) proportion')
    parser.add_argument('--path-history', default='path.history.db',
            help='input path history as an sqlite3 database file')
    parser.add_argument('--rates', default='rate.matrix.db',
            help=('continuous-time Markov substitution process '
                'as an sqlite3 database file'))
    parser.add_argument('--partition', default='partition.db',
            help=('a partition of the primary states '
                'as an sqlite3 database file'))
    main(parser.parse_args())

