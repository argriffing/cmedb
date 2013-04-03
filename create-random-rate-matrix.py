"""
Create a random rate matrix for a continuous-time Markov process.

Currently, only create time-reversible rate matrices.
"""

import argparse
import sqlite3

import numpy as np

import cmedbutil


def main(args):

    # record the number of states
    nstates = args.nstates

    # construct the rate matrix and stationary distribution
    B = np.random.exponential(scale=1.0, size=(nstates, nstates))
    S = B + B.T
    pre_distn = np.random.exponential(scale=1.0, size=(nstates,))
    distn = pre_distn / np.sum(pre_distn)
    pre_Q = S * distn
    Q = pre_Q - np.diag(np.sum(pre_Q, axis=1))

    # normalize the expected rate if requested
    if args.expected_rate:
        rates = -np.diag(Q)
        observed_rate = np.dot(rates, distn)
        Q *= (args.expected_rate / observed_rate)

    # check time-reversible rate matrix invariants
    cmedbutil.assert_stochastic_vector(distn)
    cmedbutil.assert_rate_matrix(Q)
    cmedbutil.assert_equilibrium(Q, distn)
    cmedbutil.assert_detailed_balance(Q, distn)

    # create or open the database for output
    conn = sqlite3.connect('random.rate.matrix.db')
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

    # populate the rate matrix table
    for i in range(nstates):
        for j in range(nstates):
            if i != j:
                cursor.execute(
                        'insert into rates values (?, ?, ?)',
                        (i, j, Q[i, j]),
                        )
    conn.commit()

    # populate the states table
    for i in range(nstates):
        cursor.execute('insert into states values (?, ?)', (i, str(i)))
    conn.commit()

    # populate the state distribution table
    for pair in enumerate(distn):
        cursor.execute('insert into distn values (?, ?)', pair)
    conn.commit()

    # close the database connection
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nstates', type=cmedbutil.pos_int, default=4,
            help='use this many states')
    #parser.add_argument('--reversible', action='store_true',
            #help=('require the corresponding Markov process '
                #'to be time-reversible'))
    parser.add_argument('--expected-rate', default=1.0,
            type=cmedbutil.pos_float,
            help='expected substitutions per time unit')
    main(parser.parse_args())

