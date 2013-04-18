"""
Use dynamic programming to compute a complicated likelihood.

The likelihood is over a sampled state history on an alignment of leaf states.
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
import mmpp


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
        part, partition, distn, dg, G_dag,
        rate_on, rate_off):
    """
    This uses more-clever-than-brute force likelihood calculation.
    In particular it uses dynamic programming or memoization or whatever.
    @param part: the part of the partition defining the current blink thred
    @param partition: a map from primary state to part
    @param distn: map from primary state to equilibrium probability
    @param dg: sparse primary state rate matrix as weighted directed networkx
    @param G_dag: directed phylogenetic tree with blen and state edge values
    @param rate_on: a blink rate
    @param rate_off: a blink rate
    @return: log likelihood
    """

    # Beginning at the leaves and working toward the root,
    # compute subtree likelihoods conditional on each blink state.
    v_to_b_to_lk = defaultdict(dict)

    # Initialize the likelihood map of each leaf vertex.
    leaf_set = set(v for v in G_dag if G_dag.degree(v) == 1)
    for leaf in leaf_set:

        # These likelihoods are allowed to be 1 even when
        # the leaf blink state conflicts with the primary state at the leaf,
        # because the conflicts will be handled at the edge level
        # rather than at the vertex level.
        v_to_b_to_lk[leaf][0] = 1.0
        v_to_b_to_lk[leaf][1] = 1.0

    # Work towards the root.
    for v in reversed(nx.topological_sort(G_dag)):
        if v in leaf_set:
            continue

        # prepare to multiply by likelihoods for each successor branch
        v_to_b_to_lk[v][0] = 1.0
        v_to_b_to_lk[v][1] = 1.0
        for succ in G_dag.successors(v):

            # Get the primary state of this segment,
            # and get its corresponding partition part.
            pri_state = G_dag[v][succ]['state']
            blen = G_dag[v][succ]['blen']
            pri_part = partition[part]

            # Get the conditional rate of turning off the blinking.
            # This is zero if the primary state corresponds to the
            # blink thread state, and otherwise it is rate_off.
            if partition[pri_state] == part:
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
            for sink in dg.successors(pri_state):
                rate = dg[pri_state][sink]['weight']
                if partition[sink] == part:
                    rate_absorb += rate

            # Construct the micro rate matrix and transition matrix.
            P_micro = mmpp.get_mmpp_block(
                    conditional_rate_on,
                    conditional_rate_off,
                    rate_absorb,
                    blen,
                    )
            """
            Q_micro_slow = get_micro_rate_matrix(
                    conditional_rate_off, conditional_rate_on, rate_absorb)
            P_micro_slow = scipy.linalg.expm(Q_micro_slow * blen)[:2, :2]
            if not np.allclose(P_micro, P_micro_slow):
                raise Exception((P_micro, P_micro_slow))
            """

            # Get the likelihood using the v_to_b_to_lk map.
            lk_branch = {}
            lk_branch[0] = 0.0
            lk_branch[1] = 0.0
            for ba, bb in product((0, 1), repeat=2):
                lk_transition = 1.0
                if partition[pri_state] == part:
                    if not (ba and bb):
                        lk_transition *= 0.0
                lk_transition *= P_micro[ba, bb]
                lk_rest = v_to_b_to_lk[succ][bb]
                lk_branch[ba] += lk_transition * lk_rest

            # Multiply by the likelihood associated with this branch.
            v_to_b_to_lk[v][0] *= lk_branch[0]
            v_to_b_to_lk[v][1] *= lk_branch[1]

    # get the previously arbitrarily chosen phylogenetic root vertex
    root = nx.topological_sort(G_dag)[0]

    # get the primary state at the root
    root_successor = G_dag.successors(root)[0]
    initial_primary_state = G_dag[root][root_successor]['state']

    # The initial distribution contributes to the likelihood.
    if partition[initial_primary_state] == part:
        initial_proportion_off = 0.0
        initial_proportion_on = 1.0
    else:
        initial_proportion_off = rate_off / float(rate_on + rate_off)
        initial_proportion_on = rate_on / float(rate_on + rate_off)
    path_likelihood = 0.0
    path_likelihood += initial_proportion_off * v_to_b_to_lk[root][0]
    path_likelihood += initial_proportion_on * v_to_b_to_lk[root][1]

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


def pick_arbitrary_root(G):
    """
    The vertex to be picked as the root should meet some criteria.
    The first criterion is that its degree is at least 2.
    The second criterion is that all edges adjacent to the root
    vertex should share the same primary process state.
    @param G: undirected networkx graph with state annotation of edges
    @return: an arbitrary vertex qualified to act as the root
    """
    for v in G:
        neighbors = G.neighbors(v)
        if len(neighbors) > 1:
            neighbor_edge_states = [G[v][n]['state'] for n in neighbors]
            if len(set(neighbor_edge_states)) == 1:
                return v
    raise Exception(
            'failed to find any vertex '
            'qualified to act as the root')


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
            if args.verbose:
                print history, offset
            # initialize the networkx undirected graph
            G = nx.Graph()
            for history, offset, segment, va, vb, blen, state in offset_group:
                # add an edge to the undirected graph
                G.add_edge(va, vb, blen=blen, state=state)
            # construct a directed acyclic graph with the arbitrary root
            G_dag = nx.bfs_tree(G, pick_arbitrary_root(G))
            for a, b in G_dag.edges():
                G_dag[a][b]['blen'] = G[a][b]['blen']
                G_dag[a][b]['state'] = G[a][b]['state']
            # add the log likelihood contribution of the primary thread
            ll = get_primary_log_likelihood(distn, dg, G_dag)
            log_likelihood += ll
            if args.verbose:
                print 'll contribution of primary process:', ll
            # add the log likelihood contribution of the blink threads
            ll = 0.0
            for part in range(nparts):
                ll += get_dynamic_blink_thread_log_likelihood(
                        part, partition, distn, dg, G_dag,
                        args.rate_on, args.rate_off)
            if args.verbose:
                print 'll contribution of blink process:', ll
                print
            log_likelihood += ll
        history_ll_pair = (history, log_likelihood)
        if args.verbose:
            print history_ll_pair
            print
        history_ll_pairs.append(history_ll_pair)

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
    parser.add_argument('-v', '--verbose', action='store_true',
            help='spam more text')
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

