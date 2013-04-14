"""
Use dynamic programming to compute a complicated likelihood.

The likelihood is over a sampled state history on a tree.
Not only are the states at the leaves assumed to be known,
but so are the states at the internal vertices
and also at all points along the branches.
But there is also an unobserved continuous process
consisting of the tolerated/untolerated status of each amino acid.
Because there are 20 amino acids,
the size of the state space of this unobserved process is n = 2^20
(at any time, each tolerance can be on or off),
which is too large for operations that cost n^3 such as matrix exponential.
But because of the structure of the unobserved process,
it is possible to break the likelihood integration over
the unobserved observed amino acid tolerance states
into 20 separate dynamic programming calculations over the tree.
.
A more advanced likelihood implementation will allow
more data than just a sampled history on the tree.
It would also allow some of the amino tolerances at the leaves
of the tree to be observed.
But this is not yet implemented.
.
Other unimplemented extensions to this blinking process model
include an unobserved synonymous substitution tolerance state,
and also amino acid dependent rates of blinking on and off.
.
The implementation could possibly be sped up by
using an explicit solution of the expm of the micro rate matrix,
possibly implemented in cython.
Or alternatively this whole script could be rewritten in C.
"""

import argparse
import sqlite3
import math
from itertools import permutations
from itertools import product
from itertools import groupby
from collections import defaultdict
import functools

import numpy as np
import networkx as nx
import scipy.linalg

import cmedbutil


# XXX most of this script is copypasted;
# the reusable functions should be moved into a separate module.


#XXX copypasted
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


def get_dynamic_blink_thread_log_likelihood(
        part, partition, distn, dg, path_history,
        rate_on, rate_off):
    """
    This uses more-clever-than-brute force likelihood calculation.
    In particular it uses dynamic programming or memoization or whatever.
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

    # At each interesting point,
    # compute the likelihood of the remaining path for (conditional on)
    # each of the two possible current blink states.
    lk_table = np.zeros((npoints, 2), dtype=float)
    lk_table[-1] = 1
    for segment_index in reversed(range(nsegments)):
        segment, primary_state, duration = path_history[segment_index]

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

        # Get the likelihood using the table.
        for ba, bb in product((0, 1), repeat=2):
            lk_transition = 1.0
            if partition[primary_state] == part:
                if not (ba and bb):
                    lk_transition *= 0.0
            lk_transition *= P_micro[ba, bb]
            lk_rest = lk_table[segment_index+1, bb]
            lk_table[segment_index, ba] += lk_transition * lk_rest

    # slowly and inefficiently get the first primary state
    seg_seq, primary_state_seq, duration_seq = zip(*path_history)
    initial_primary_state = primary_state_seq[0]

    # The initial distribution contributes to the likelihood.
    if partition[initial_primary_state] == part:
        initial_proportion_off = 0.0
        initial_proportion_on = 1.0
    else:
        initial_proportion_off = rate_off / float(rate_on + rate_off)
        initial_proportion_on = rate_on / float(rate_on + rate_off)
    path_likelihood = 0.0
    path_likelihood += initial_proportion_off * lk_table[0, 0]
    path_likelihood += initial_proportion_on * lk_table[0, 1]

    # Report the log likelihood.
    return math.log(path_likelihood)


def get_primary_log_likelihood(distn, dg, G_dag):
    """
    This includes the initial state and the transitions.
    It does not include the dwell times,
    because these contributions are added by the blinking threads.
    @param distn: map from primary state to equilibrium probability
    @param dg: sparse primary state rate matrix as weighted directed networkx
    @param G_dag: directed phylogenetic tree with blen and state edge values
    @return: log likelihood
    """

    # get the previously arbitrarily chosen phylogenetic root vertex
    root = nx.topological_sort(G_dag)[0]

    # initialize log likelihood
    log_likelihood = 0.0

    # get the primary state at the root
    root_successor = G_dag.successors(root)[0]
    initial_primary_state = G_dag[root][root_successor]['state']

    # add the contribution of the equilibrium distribution
    log_likelihood += math.log(distn[initial_primary_state])

    # add the contribution of the primary state transitions
    for v in G_dag:
        preds = G_dag.predecessors(v)
        succs = G_dag.successors(v)
        if len(preds) == 1 and len(succs) == 1:
            pred = preds[0]
            succ = succs[0]
            a = G_dag[pred][v]['state']
            b = G_dag[v][succ]['state']

            # Note that the states can be the same
            # if there was a degree two vertex in the original tree,
            # for example if the original tree was rooted.
            if a != b:
                rate = dg[a][b]['weight']
                log_likelihood += math.log(rate)

    # return the log likelihood
    return log_likelihood


#XXX massive copypasting from a similarly named script
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

    # Open the primary state tree history from the sqlite3 database file.
    conn = sqlite3.connect(args.histories)
    cursor = conn.cursor()

    # Read rows of sampled tree history from the database,
    # being careful to not read too much into RAM at once.
    # Build the list of log likelihoods.
    # There should be one log likelihood for each history.
    cursor.execute(
            'select history, offset, segment, va, vb, blen, state '
            'from histories '
            'order by history, offset')
    history_ll_pairs = []
    first_element = cmedbutil.first_element
    second_element = functools.partial(cmedbutil.nth_element, 1)
    for history, history_group in groupby(cursor, first_element):
        log_likelihood = 0.0
        for offset, offset_group in groupby(history_group, second_element):
            # initialize the networkx undirected graph
            G = nx.Graph()
            for history, offset, segment, va, vb, blen, state in offset_group:
                # add an edge to the undirected graph
                G.add_edge(va, vb, blen=blen, state=state)
            # choose a high-degree internal vertex as the root
            deg_v_pairs = [(v, G.degree(v)) for v in G]
            root_degree, root = max(deg_v_pairs)
            # construct a directed acyclic graph with the arbitrary root
            G_dag = nx.bfs_tree(G, root)
            for a, b in G_dag.edges():
                G_dag[a][b]['blen'] = G[a][b]['blen']
                G_dag[a][b]['state'] = G[a][b]['state']
            # add the log likelihood contribution of the primary thread
            log_likelihood += get_primary_log_likelihood(
                    distn, dg, G_dag)
            # add the log likelihood contribution of the blink threads
            #XXX
            """
            for part in range(nparts):
                log_likelihood += get_dynamic_blink_thread_log_likelihood(
                        part, partition, distn, dg, G_dag,
                        args.rate_on, args.rate_off)
            """
        history_ll_pairs.append((history, log_likelihood))

    # close the input database of tree histories
    conn.close()

    # open the database file for the output log likelihoods
    conn = sqlite3.connect(args.outfile)
    cursor = conn.cursor()

    # initialize the log likelihoods table
    s = (
            'create table if not exists log_likelihoods ('
            'history integer, '
            'log_likelihood real, '
            'primary key (history))')
    cursor.execute(s)
    conn.commit()

    # write the log likelihoods into the table
    for t in history_ll_pairs:
        s = 'insert into log_likelihoods values (?, ?)'
        cursor.execute(s, t)
    conn.commit()

    # close the database file for the patterns
    conn.close()


if __name__ == '__main__':

    # begin to define the command line arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rates', default='rate.matrix.db',
            help='time-reversible rate matrix as an sqlite3 database file')
    parser.add_argument('--histories', default='tree.histories.db',
            help='input tree histories as an sqlite3 database file')
    parser.add_argument('--outfile', default='log.likelihoods.db',
            help='output log likelihoods in sqlite3 format')

    # these args are specific to the blinking process
    parser.add_argument('--rate-on',
            type=cmedbutil.nonneg_float,
            help='rate at which blink states change from off to on')
    parser.add_argument('--rate-off',
            type=cmedbutil.nonneg_float,
            help='rate at which blink states change from on to off')
    parser.add_argument('--partition', default='partition.db',
            help=('a partition of the primary states '
                'as an sqlite3 database file'))

    # call the main function
    main(parser.parse_args())

