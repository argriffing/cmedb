"""
Construct a fast blinking process limit of an existing rate matrix.

The set of states is assumed to be partitioned.
The modification of the rate matrix
involves only reducing some rates and adjusting the distribution.
"""

import argparse
import sqlite3
from itertools import permutations

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



def main(args):

    # Read the sparse rate matrix from a database file.
    # Also read the state, name pairs for copypasting to the output matrix.
    conn = sqlite3.connect(args.rates)
    cursor = conn.cursor()
    distn, dg = get_sparse_rate_matrix_info(cursor)
    state_name_pairs = list(cursor.execute('select state, name from states'))
    conn.close()

    # Read the partition from a database file.
    # The hidden blinking process controls the state transition
    # of the primary process according to this partition
    # of primary process states.
    conn = sqlite3.connect(args.partition)
    cursor = conn.cursor()
    partition = dict(cursor.execute('select state, part from partition'))
    conn.close()

    # adjust the rates
    for a in dg:
        for b in dg.successors(a):
            rate = dg[a][b]['weight']
            if partition[a] == partition[b]:
                rate *= args.syn_proportion_on
            else:
                rate *= args.proportion_on
            dg[a][b]['weight'] = rate

    # the stationary distribution does not change

    # create or open the rate matrix database for writing
    conn = sqlite3.connect(args.outfile)
    cursor = conn.cursor()

    # create the tables
    cursor.execute(
            'create table rates ('
            'source integer, sink integer, rate real, '
            'primary key (source, sink))')
    cursor.execute(
            'create table states ('
            'state integer, name text, '
            'primary key (state))')
    cursor.execute(
            'create table distn ('
            'state integer, prob real, '
            'primary key (state))')
    conn.commit()

    # populate the output rate matrix table
    for a in dg:
        for b in dg.successors(a):
            rate = dg[a][b]['weight']
            if rate:
                s = 'insert into rates values (?, ?, ?)'
                t = (a, b, rate)
                cursor.execute(s, t)
    conn.commit()

    # populate the states table
    for pair in state_name_pairs:
        cursor.execute('insert into states values (?, ?)', pair)
    conn.commit()

    # populate the state distribution table
    for item in distn.items():
        cursor.execute('insert into distn values (?, ?)', item)
    conn.commit()

    # close the database connection
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rates', default='rate.matrix.db',
            help=('input continuous-time Markov chain rate matrix '
                'as an sqlite3 database file'))
    parser.add_argument('--partition', default='partition.db',
            help=('input partition of the primary states '
                'as an sqlite3 database file'))
    parser.add_argument('--proportion-on', type=cmedbutil.nonneg_float,
            help='part toleration probability')
    parser.add_argument('--syn-proportion-on', type=cmedbutil.nonneg_float,
            help='synonymous substitution toleration probability')
    parser.add_argument('--outfile', default='fast.blink.rate.matrix.db',
            help='output rate matrix in sqlite3 format')
    main(parser.parse_args())

