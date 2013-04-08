"""
Compute a blinking-process probability of primary state endpoints.

This is on a path, not a tree.
Only the initial and final primary states are observed.
The primary state is not observed along the path,
and the blinking process is not observed anywhere.
This script is slow, and it computes an exact probability
by summing over the uninknown blinking process initial state distribution
and integrating over all of the unknown blinking process transitions
and over all of the primary path transitions.
This is not feasible for actual genetic codes.
.
The practical purpose of this script is to work in tandem with
artificial genetic codes for the purposes of testing "Monte Carlo
likelihood ratio" (MCLR) implementations.
If Monte Carlo likelihood ratio estimation seems to converge to the
exact likelihood ratios for small artificial genetic codes,
then the MCLR implementation would presumably not be buggy and
could be used for inference with non-toy genetic codes.
.
Test using a toy rate matrix with random partitions.
Possibly turn this into a makefile for testing purposes.
$ python create-random-rate-matrix.py
--nstates=8 --outfile=toy.rate.matrix.db
$ python bp-random-state-partition.py
--outfile=toy.partition.db --nparts=4 --rates=toy.rate.matrix.db
$ python ctmc-segment-bridge-sampling.py
--initial=0 --final=7 --method=modified-rejection --elapsed=5
--rates=toy.rate.matrix.db --nsamples=100 --table=histories
--outfile=toy.histories.db
$ python path-histories-blink-likelihoods.py
--rates=toy.rate.matrix.db --histories=toy.histories.db
--rate-on=1 --rate-off=3 --partition=toy.partition.db
--outfile=toy.blink.log.likelihoods.db
$ python path-histories-likelihoods.py
--rates=toy.rate.matrix.db --histories=toy.histories.db
--outfile=toy.reference.log.likelihoods.db
$ python monte-carlo-likelihood-ratio.py
--numerator-log-likelihood=toy.blink.log.likelihoods.db
--denominator-log-likelihood=toy.reference.log.likelihoods.db 
$ python ctmc-endpoint-probability.py
--elapsed=5 --initial=0 --final=7 --rates=toy.rate.matrix.db
$ python bp-endpoint-probability.py
--elapsed=5 --initial=0 --final=7 --rates=toy.rate.matrix.db
--partition=toy.partition.db --rate-on=1 --rate-off=3
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


def create_compound_rate_matrix(rate_off, rate_on, partition, distn, dg):
    """
    Construct the augmented rate matrix.
    The augmented rate matrix includes states corresponding
    to all combinations of primary and primary-tolerance states.
    The output includes
    a map from compound state to primary state as a list,
    a map from compound state to blinking state vector as a list,
    a dense compound rate matrix.
    @param rate_off: rate from the 'on' state to the 'off' state
    @param rate_on: rate from the 'off' state to the 'on' state
    @param partition: a map from primary state to part
    @param distn: map from primary state to equilibrium probability
    @param dg: sparse primary state rate matrix as weighted directed networkx
    @return: primary map, blink map, compound rate matrix
    """

    # get the set of indices of parts
    parts = set(partition.values())
    nparts = len(parts)

    # define the compound states
    compound_to_primary = []
    compound_to_blink = []
    for primary in distn:
        for blink_array in product((0, 1), repeat=nparts):
            compound_to_primary.append(primary)
            compound_to_blink.append(blink_array)

    # count the compound states
    nprimary = len(distn)
    ncompound = nprimary * (2 ** nparts)

    # define the off-diagonal entries of the dense compound rate matrix
    pre_Q = np.zeros((ncompound, ncompound), dtype=float)
    for i in range(ncompound):
        pri_i = compound_to_primary[i]
        blink_i = compound_to_blink[i]
        for j in range(ncompound):

            # we are only concerned with compound state changes
            if i == j:
                continue

            # decompose the compound states
            pri_j = compound_to_primary[j]
            blink_j = compound_to_blink[j]

            # precompute the number of blink transitions
            blink_hdist = hamming_distance(blink_i, blink_j)

            # we do not allow multiple simultaneous blink transitions
            if blink_hdist > 1:
                continue

            # simultaneous primary and blink state changes are not allowed
            if pri_i != pri_j and blink_hdist == 1:
                continue

            # get the primary state partitions
            part_i = partition[pri_i]
            part_j = partition[pri_j]

            # check for primary state transition
            if pri_i != pri_j:

                # primary state cannot change unless both blink states are on
                if not all((
                        blink_i[part_i],
                        blink_i[part_j],
                        blink_j[part_i],
                        blink_j[part_j],
                        )):
                    continue

                # set the compound state transition rate
                if pri_j in dg[pri_i]:
                    rate = dg[pri_i][pri_j]['weight']
                    pre_Q[i, j] = rate

            # check for single blink state transition
            if blink_hdist == 1:

                # get the blink difference and part
                diff = sum(b-a for a, b in zip(blink_i, blink_j))
                blink_part = index_of_first_difference(blink_i, blink_j)

                # set the compound state transition rate
                if diff == 1:
                    rate = rate_on
                elif diff == -1:
                    if blink_part == part_i:
                        rate = 0
                    else:
                        rate = rate_off
                else:
                    raise Exception
                pre_Q[i, j] = rate

    # define the diagonal entries of the compound rate matrix
    Q = pre_Q - np.diag(np.sum(pre_Q, axis=1))

    # return compound state transition rate information
    return compound_to_primary, compound_to_blink, Q


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
    ncompound = nprimary * (2 ** nparts)

    # get the compound state transition rate info
    compound_to_primary, compound_to_blink, Q = create_compound_rate_matrix(
            args.rate_off, args.rate_on, partition, distn, dg)

    # compute compound equilibrium distribution
    compound_distn = np.zeros(ncompound, dtype=float)
    blink_distn = np.array([
        args.rate_off / (args.rate_on + args.rate_off),
        args.rate_on / (args.rate_on + args.rate_off),
        ], dtype=float)
    for i in range(ncompound):
        pri = compound_to_primary[i]
        part = partition[pri]
        blink = compound_to_blink[i]
        prob = distn[pri]
        for blink_part, blink_state in enumerate(blink):
            if part == blink_part:
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
    for i in range(ncompound):
        pri_i = compound_to_primary[i]
        for j in range(ncompound):
            pri_j = compound_to_primary[j]
            if pri_i == args.initial and pri_j == args.final:
                joint_prob += compound_distn[i] * P[i, j]

    # report the joint probability
    print 'joint log probability:', math.log(joint_prob)
    print 'joint probability:', joint_prob



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rate-on',
            type=cmedbutil.nonneg_float,
            help='rate at which blink states change from off to on')
    parser.add_argument('--rate-off',
            type=cmedbutil.nonneg_float,
            help='rate at which blink states change from on to off')
    parser.add_argument('--elapsed',
            type=cmedbutil.nonneg_float,
            help='elapsed time between endpoints')
    parser.add_argument('--initial', type=cmedbutil.nonneg_int, required=True,
            help='initial state')
    parser.add_argument('--final', type=cmedbutil.nonneg_int, required=True,
            help='final state')
    parser.add_argument('--rates', default='rate.matrix.db',
            help=('continuous-time Markov substitution process '
                'as an sqlite3 database file'))
    parser.add_argument('--partition', default='partition.db',
            help=('a partition of the primary states '
                'as an sqlite3 database file'))
    main(parser.parse_args())

