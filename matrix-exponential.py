"""
Compute the transition matrix given the rate matrix.

This function does not currently do anything clever
with the sparsity of a rate matrix.
"""

import argparse
import sqlite3
import math
import itertools

import numpy as np
import scipy.linalg

import cmedbutil


def get_rate_matrix_info(cursor):
    """
    @param cursor: sqlite3 database cursor
    @return: sorted state list, stationary distribution, dense rate matrix
    """

    # get the sorted list of all states
    cursor.execute(
            'select state from distn '
            'union '
            'select source from rates '
            'union '
            'select sink from rates '
            'union '
            'select state from states '
            )
    states = sorted(t[0] for t in cursor)

    # count the states
    nstates = len(states)

    # define the map from state to rate matrix index
    s_to_i = dict((s, i) for i, s in enumerate(states))

    # construct the rate matrix
    cursor.execute('select source, sink, rate from rates')
    pre_Q = np.zeros((nstates, nstates), dtype=float)
    for si, sj, rate in cursor:
        pre_Q[s_to_i[si], s_to_i[sj]] = rate
    Q = pre_Q - np.diag(np.sum(pre_Q, axis=1))

    # construct the distribution of states at the root
    cursor.execute('select state, prob from distn')
    state_prob_pairs = list(cursor)
    distn = np.zeros(nstates)
    for state, prob in state_prob_pairs:
        distn[s_to_i[state]] = prob

    # assert that the distribution has the right form
    cmedbutil.assert_stochastic_vector(distn)

    # assert that the rate matrix is actually a rate matrix
    cmedbutil.assert_rate_matrix(Q)

    # assert that the distribution is at equilibrium w.r.t. the rate matrix
    cmedbutil.assert_equilibrium(Q, distn)

    # assert that the detailed balance equations are met
    cmedbutil.assert_detailed_balance(Q, distn)

    # return the validated inputs describing the stochastic process
    return states, distn, Q


def main(args):

    # read and validate the rate matrix info from the sqlite3 database file
    conn = sqlite3.connect(args.rates)
    cursor = conn.cursor()
    states, distn, Q = get_rate_matrix_info(cursor)
    conn.close()

    # extract the amount of time along the path
    T = args.elapsed

    # Compute the matrix exponential M.
    # This matrix is often represented as P but the Holmes-Rubin (2002)
    # paper uses the notation M.
    M = scipy.linalg.expm(Q*T)

    # create the output database file and initialize the cursor
    conn = sqlite3.connect(args.outfile)
    cursor = conn.cursor()

    # define the transition matrix table
    s = (
            'create table if not exists transitions ('
            'source integer, '
            'sink integer, '
            'prob real, '
            'primary key (source, sink))')
    cursor.execute(s)
    conn.commit()

    # populate the transitions table
    for i, source_state in enumerate(states):
        for j, sink_state in enumerate(states):
            s = 'insert into transitions values (?, ?, ?)'
            t = (source_state, sink_state, M[i, j])
            cursor.execute(s, t)
    conn.commit()

    # close the output database connection
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--elapsed', type=cmedbutil.pos_float, default=1.0,
            help='elapsed time')
    parser.add_argument('--rates', default='rate.matrix.db',
            help='time-reversible rate matrix as an sqlite3 database file')
    parser.add_argument('--outfile', default='transition.matrix.db',
            help='output expectations as an sqlite3 database file')
    main(parser.parse_args())

