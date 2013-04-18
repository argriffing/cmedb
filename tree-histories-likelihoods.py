"""
Compute a continuous time Markov chain likelihood for histories on trees.

Strangely this script is being written after the more complicated
blinking process likelihood script.
This script does not require dynamic programming;
it just requires counting things and computing scalar exponentials.
It does not require matrix exponentials.
However it does require the rate matrix to be time-reversible.
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


# XXX directly copypasted
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

    # read the sparse rate matrix from a database file
    conn = sqlite3.connect(args.rates)
    cursor = conn.cursor()
    distn, dg = get_sparse_rate_matrix_info(cursor)
    conn.close()

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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--verbose', action='store_true',
            help='spam more text')
    parser.add_argument('--rates', default='rate.matrix.db',
            help='time-reversible rate matrix as an sqlite3 database file')
    parser.add_argument('--histories', default='tree.histories.db',
            help='input tree histories as an sqlite3 database file')
    parser.add_argument('--outfile', default='log.likelihoods.db',
            help='output log likelihoods in sqlite3 format')
    main(parser.parse_args())

