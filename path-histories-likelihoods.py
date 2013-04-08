"""
Compute likelihoods for multiple sampled path histories.

This script should not compute the expm directly.
"""

import argparse
import sqlite3
import math
from itertools import permutations
from collections import defaultdict

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


#XXX massive copypasting from a similarly named script
def main(args):

    # read the sparse rate matrix from a database file
    conn = sqlite3.connect(args.rates)
    cursor = conn.cursor()
    distn, dg = get_sparse_rate_matrix_info(cursor)
    conn.close()

    # read the primary state path history from the sqlite3 database file
    conn = sqlite3.connect(args.histories)
    cursor = conn.cursor()
    cursor.execute(
            'select history, segment, state, blen '
            'from histories '
            'order by history, segment')
    data = list(cursor)
    conn.close()

    # put the paths into a dictionary
    history_to_path = defaultdict(list)
    for history, segment, state, blen in data:
        path_history = (segment, state, blen)
        history_to_path[history].append(path_history)

    # compute the path history log likelihoods
    history_ll_pairs = []
    for history, path_history in history_to_path.items():

        # initialize
        log_likelihood = 0.0
        seg_seq, state_seq, duration_seq = zip(*path_history)

        # add the contribution of the equilibrium distribution
        initial_state = state_seq[0]
        log_likelihood += math.log(distn[initial_state])

        # add the contribution of the primary state transitions
        for a, b in cmedbutil.pairwise(state_seq):
            rate = dg[a][b]['weight']
            log_likelihood += math.log(rate)

        # add the contribution of the dwell times
        for segment, state, duration in path_history:
            total_rate = sum(
                    dg[state][b]['weight'] for b in dg.successors(state))
            log_likelihood -= duration * total_rate

        # add the (history, log_likelihood) to the list
        history_ll_pairs.append((history, log_likelihood))

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
    parser.add_argument('--rates', default='rate.matrix.db',
            help='time-reversible rate matrix as an sqlite3 database file')
    parser.add_argument('--histories', default='path.histories.db',
            help='input path histories as an sqlite3 database file')
    parser.add_argument('--outfile', default='log.likelihoods.db',
            help='output log likelihoods in sqlite3 format')
    main(parser.parse_args())

