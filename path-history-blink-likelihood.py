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
from itertools import product

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


def get_partially_observed_blink_thread_log_likelihood(
        part, partition, distn, dg, path_history,
        rate_on, rate_off,
        endpoint_assignment,
        ):
    """
    @param part: the part of the partition defining the current blink thread
    @param partition: a map from primary state to part
    @param distn: map from primary state to equilibrium probability
    @param dg: sparse primary state rate matrix as weighted directed networkx
    @param path_history: triples of (segment, state, blen)
    @param rate_on: a blink rate
    @param rate_off: a blink rate
    @param endpoint_assignment: an assignment of blink states at primary trans
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
        if not endpoint_assignment[0]:
            return -float('Inf')
    else:
        blink_distn = np.array([
            rate_off / float(rate_on + rate_off),
            rate_on / float(rate_on + rate_off),
            ], dtype=float)
        log_blink_distn = np.log(blink_distn)
        finite_path_ll += log_blink_distn[endpoint_assignment[0]]

    # multiply the likelihood across all segments along the thread
    for i in range(nsegments):

        # extract the primary and blink state
        segment, primary_state, duration = path_history[i]
        ba = endpoint_assignment[i]
        bb = endpoint_assignment[i+1]

        # If the blinking state is 'off' at an endpoint
        # adjacent to a segment with a compatible primary state
        # then set the likelihood to zero.
        if partition[primary_state] == part:
            if not (ba and bb):
                return -float('Inf')

        # Get the conditional rate of turning off the blinking.
        # This is zero if the primary state corresponds to the
        # blink thread state, and otherwise it is rate_off.
        if partition[primary_state] == part:
            conditional_rate_off = 0.0
        else:
            conditional_rate_off = rate_off

        # Get the conditional rate of turning on the blinking.
        # This is always rate_on.
        conditional_rate_on = rate_on

        # Get the absorption rate.
        # This is the sum of primary transition rates
        # into the part that corresponds to the current blink thread state.
        rate_absorb = 0.0
        for sink in dg.successors(primary_state):
            rate = dg[primary_state][sink]['weight']
            if partition[sink] == part:
                rate_absorb += rate

        # Construct the micro rate matrix and transition matrix.
        Q_micro = get_micro_rate_matrix(
                conditional_rate_off, conditional_rate_on, rate_absorb)
        P_micro = scipy.linalg.expm(Q_micro * duration)
        P_micro_log = np.log(P_micro)

        # Contribute to the likelihood.
        finite_path_ll += P_micro_log[ba, bb]

    # Return the finite path likelihood for this blinking thread.
    return finite_path_ll


def get_blink_thread_log_likelihood(
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
    @return: log likelihood
    """

    # count the number of segments and the number of segment endpoints
    nsegments = len(path_history)
    npoints = nsegments + 1

    # initialize the likelihood
    log_likelihoods = []

    # sum the likelihood over all possible states at the segment endpoints
    for endpoint_assignment in product((0, 1), repeat=npoints):
        log_likelihood = get_partially_observed_blink_thread_log_likelihood(
                part, partition, distn, dg, path_history,
                rate_on, rate_off,
                endpoint_assignment)
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

    # check blink rates
    if (args.rate_on + args.rate_off) in (0, float('Inf')):
        raise NotImplementedError(
                'for extreme limits of blinking rates, '
                'use a different script to compute the likelihood')

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
                args.rate_on, args.rate_off)

    # report the likelihood
    print 'log likelihood:', log_likelihood
    print 'likelihood:', math.exp(log_likelihood)



if __name__ == '__main__':

    # define some methods
    method_choices = (
            'brute',
            'dynamic',
            )

    # define the command line parameters
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rate-on',
            type=cmedbutil.nonneg_float,
            help='rate at which blink states change from off to on')
    parser.add_argument('--rate-off',
            type=cmedbutil.nonneg_float,
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

