"""
Compute a blinking-process probability of primary state endpoints.

This script has been copypasted and is expected to have been modified
to allow a blinking state that enables/disables synonymous substitutions.
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



def hamming_distance(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)

def index_of_first_difference(a, b):
    for i, (x, y) in enumerate(zip(a, b)):
        if x != y:
            return i


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


def create_compound_rate_matrix(
        rate_off, rate_on,
        syn_rate_off, syn_rate_on,
        partition, distn, dg):
    """
    Construct the augmented rate matrix.
    The augmented rate matrix includes states corresponding
    to all combinations of primary, syn-tolerance, primary-tolerance states.
    The output includes
    a compound state list
    a dense compound rate matrix.
    @param rate_off: rate from the 'on' state to the 'off' state
    @param rate_on: rate from the 'off' state to the 'on' state
    @param syn_rate_off: rate from the syn 'on' state to the 'off' state
    @param syn_rate_on: rate from the syn 'off' state to the 'on' state
    @param partition: a map from primary state to part
    @param distn: map from primary state to equilibrium probability
    @param dg: sparse primary state rate matrix as weighted directed networkx
    @return: primary map, blink map, compound rate matrix
    """

    # get the set of indices of parts
    parts = set(partition.values())
    nparts = len(parts)

    # Define the ordered list of hashable compound states.
    # Because the compound states are hashable,
    # the inverse of the sequence can be stored as a python dict.
    compound_states = []
    for primary in distn:
        for syn_blink in (0, 1):
            for blink_array in product((0, 1), repeat=nparts):
                compound_state = (primary, syn_blink, blink_array)
                compound_states.append(compound_state)
    compound_to_i = dict((s, i) for i, s in enumerate(compound_states))
    
    # count the compound states
    nprimary = len(distn)
    ncompound = len(compound_states)

    # define the off-diagonal entries of the dense compound rate matrix
    pre_Q = np.zeros((ncompound, ncompound), dtype=float)
    for i, (pri_i, syn_i, blink_i) in enumerate(compound_states):

        # define the partition part of the primary source state
        part_i = partition[pri_i]

        # Explicitly construct the sink states
        # that can be directly reached from the source state.
        # First add rates corresponding to allowed primary state transitions.
        # Transitions to blinked-off parts are not allowed.
        # Synonymous transitions are not allowed if the syn blink state is off.
        for pri_j in dg.successors(pri_i):
            part_j = partition[pri_j]
            if not blink_i[part_j]:
                continue
            if part_i == part_j and not syn_i:
                continue
            sink = (pri_j, syn_i, blink_i)
            j = compound_to_i[sink]
            pre_Q[i, j] = dg[pri_i][pri_j]['weight']

        # Add rates corresponding to synonymous blink transitions.
        if syn_i:
            j = compound_to_i[pri_i, 0, blink_i]
            pre_Q[i, j] = syn_rate_off
        else:
            j = compound_to_i[pri_i, 1, blink_i]
            pre_Q[i, j] = syn_rate_on

        # Add rates corresponding to allowed part blink transitions.
        # Blinking off the current part tolerance is not allowed.
        for part_k, blink_k in enumerate(blink_i):

            # It is not allowed to change the tolerance state
            # corresponding to the current primary part.
            if part_i == part_k:
                continue

            # Blink the tolerance state on or off.
            blink_j = list(blink_i)
            if blink_k:
                blink_j[part_k] = 0
                rate = rate_off
            else:
                blink_j[part_k] = 1
                rate = rate_on
            j = compound_to_i[pri_i, syn_i, tuple(blink_j)]
            pre_Q[i, j] = rate

    # define the diagonal entries of the compound rate matrix
    Q = pre_Q - np.diag(np.sum(pre_Q, axis=1))

    # return compound state transition rate information
    return compound_states, Q


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

    # count states
    parts = set(partition.values())
    nparts = len(parts)
    nprimary = len(distn)

    # get the compound state transition rate info
    compound_states, Q = create_compound_rate_matrix(
            args.rate_off, args.rate_on,
            args.syn_rate_off, args.syn_rate_on,
            partition, distn, dg)
    ncompound = len(compound_states)

    # compute compound equilibrium distribution
    compound_distn = np.zeros(ncompound, dtype=float)
    blink_distn = np.array([
        args.rate_off / (args.rate_on + args.rate_off),
        args.rate_on / (args.rate_on + args.rate_off),
        ], dtype=float)
    syn_blink_distn = np.array([
        args.syn_rate_off / (args.syn_rate_on + args.syn_rate_off),
        args.syn_rate_on / (args.syn_rate_on + args.syn_rate_off),
        ], dtype=float)
    for i, (pri, syn, blink) in enumerate(compound_states):
        prob = 1.0
        prob *= distn[pri]
        prob *= syn_blink_distn[syn]
        for blink_part, blink_state in enumerate(blink):
            if partition[pri] == blink_part:
                if not blink_state:
                    prob *= 0
            else:
                prob *= blink_distn[blink_state]
        compound_distn[i] = prob

    # check the rate matrix and the distribution
    cmedbutil.assert_stochastic_vector(compound_distn)
    cmedbutil.assert_rate_matrix(Q)
    cmedbutil.assert_equilibrium(Q, compound_distn)
    cmedbutil.assert_detailed_balance(Q, compound_distn)

    # compute compound matrix exponential
    P = scipy.linalg.expm(args.elapsed * Q)

    # Compute the joint probability of the initial and final primary states.
    joint_prob = 0.0
    for i, (pri_i, syn_i, blink_i) in enumerate(compound_states):
        for j, (pri_j, syn_j, blink_j) in enumerate(compound_states):
            if pri_i == args.initial and pri_j == args.final:
                joint_prob += compound_distn[i] * P[i, j]

    # report the joint probability
    print 'joint log probability:', math.log(joint_prob)
    print 'joint probability:', joint_prob

    # if in a spammy mood, report the compound equilibrium distribution
    if args.verbose:
        print 'equilibrium distribution:'
        print compound_distn



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--verbose', action='store_true',
            help='spam more text')
    parser.add_argument('--rates', default='rate.matrix.db',
            help=('continuous-time Markov substitution process '
                'as an sqlite3 database file'))
    parser.add_argument('--partition', default='partition.db',
            help=('a partition of the primary states '
                'as an sqlite3 database file'))

    # define the path endpoints and length
    parser.add_argument('--initial', type=cmedbutil.nonneg_int, required=True,
            help='initial state')
    parser.add_argument('--final', type=cmedbutil.nonneg_int, required=True,
            help='final state')
    parser.add_argument('--elapsed',
            type=cmedbutil.nonneg_float,
            help='elapsed time between endpoints')

    # part tolerance enable/disable rates
    parser.add_argument('--rate-on',
            type=cmedbutil.nonneg_float,
            help='rate at which blink states change from off to on')
    parser.add_argument('--rate-off',
            type=cmedbutil.nonneg_float,
            help='rate at which blink states change from on to off')

    # synonymous substitution enable/disable rates
    parser.add_argument('--syn-rate-on',
            type=cmedbutil.nonneg_float,
            help='rate at which syn blink states change from off to on')
    parser.add_argument('--syn-rate-off',
            type=cmedbutil.nonneg_float,
            help='rate at which syn blink states change from on to off')

    main(parser.parse_args())

