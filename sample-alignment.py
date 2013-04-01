"""
Duplicate seq-gen because I am allergic to programs that other people write.

The continuous time Markov process is assumed to be time-reversible.
output schema
table alignment
taxon integer, offset integer, state integer
"""

import cmedbutil


#FIXME allow non-reversible rate matrices if the root of the tree is specified
#maybe use a different script for this.


#XXX multiple scripts use this function
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
    cmedb.assert_stochastic_vector(distn)

    # assert that the rate matrix is actually a rate matrix
    cmedb.assert_rate_matrix(Q)

    # assert that the distribution is at equilibrium w.r.t. the rate matrix
    cmedb.assert_equilibrium(Q, distn)

    # assert that the detailed balance equations are met
    cmedb.assert_detailed_balance(Q, distn)

    # return the validated inputs describing the stochastic process
    return states, distn, Q


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tree', required=True,
            help='unrooted tree as an sqlite3 database file')
    parser.add_argument('--rates', required=True,
            help=('continuous-time Markov substitution process '
                'as an sqlite3 database file'))
    parser.add_argument('--length', type=pos_int, default=10,
            help='sequence length')
    parser.add_argument('-o', '--outfile', default='alignment.sample.db',
            help='create this sqlite3 database file')
    main(parser.parse_args())

