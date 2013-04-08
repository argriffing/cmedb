"""
We are making the blinking process complicated by adding heterogeneous rates.
"""

import random
import argparse
import sqlite3

import numpy as np

import cmedbutil

#XXX not started


#XXX copypasted
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

    # construct the random partition with no empty parts
    nstates = len(states)
    nparts = args.nparts
    if nstates < nparts:
        raise Exception('expected at least as many states as parts')
    parts = [[] for i in range(nparts)]
    for i, state in enumerate(states[:nparts]):
        parts[i].append(state)
    for state in states[nparts:]:
        parts[random.randrange(nparts)].append(state)

    # open the database file for the patterns
    conn = sqlite3.connect(args.outfile)
    cursor = conn.cursor()

    # initialize the pattern multiplicities table
    s = (
            'create table if not exists partition ('
            'state integer, '
            'part integer, '
            'primary key (state))')
    cursor.execute(s)
    conn.commit()

    # write the partition information into the table
    for i, part in enumerate(parts):
        for state in part:
            s = 'insert into partition values (?, ?)'
            t = (state, i)
            cursor.execute(s, t)
    conn.commit()

    # close the database file for the patterns
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rates', default='rate.matrix.db',
            help=('continuous-time Markov substitution process '
                'as an sqlite3 database file'))
    parser.add_argument('--nparts', type=cmedbutil.pos_int, default=4,
            help='partition the state into this many parts')
    parser.add_argument('-o', '--outfile', default='blink.partition.db',
            help='create this sqlite3 database file')
    main(parser.parse_args())

