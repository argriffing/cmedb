"""
Sample the state history of a time interval given initial and final state.
"""

#XXX implement more path sampling methods
#XXX document the database formats
#XXX allow single vs multiple path samples

import functools
import argparse
import sqlite3

import numpy as np


#XXX this should go into a separate module
def nonneg_int(x):
    x = int(x)
    if x < 0:
        raise argparse.ArgumentTypeError(
                'value must be a non-negative integer')
    return x

#XXX this should go into a separate module
def pos_float(x):
    x = float(x)
    if x <= 0:
        raise argparse.ArgumentTypeError(
                'value must be a positive floating point number')
    return x

#XXX this is copypasted
def random_category(distn):
    """
    Sample from a categorical distribution.
    Note that this is not the same as random.choice(distn).
    Maybe a function like this will eventually appear
    in python or numpy or scipy.
    @param distn: categorical distribution as a stochastic vector
    @return: category index as a python integer
    """
    nstates = len(distn)
    np_index = np.dot(np.arange(nstates), np.random.multinomial(1, distn))
    return int(np_index)

#XXX this is copypasted
def decompose_rates(Q):
    """
    Break a rate matrix into two parts.
    The first part consists of the rates away from each state;
    this information is contained in the diagonal of the rate matrix.
    The second part consists of a transition matrix
    that defines the distribution over sink states conditional
    on an instantaneous change away from a given source state.
    Note that this function never requires expm of Q.
    Also, P preserves the sparsity pattern of Q.
    @param Q: rate matrix
    @return: rates, P
    """
    nstates = len(Q)
    rates = -np.diag(Q)
    P = np.array(Q)
    for i, rate in enumerate(rates):
        if rate:
            P[i, i] = 0
            P[i] /= rate
    return rates, P

#XXX this is copypasted
def gen_branch_history_sample(state_in, blen_in, rates, P):
    """
    Path sampling along a branch with a known initial state.
    Yield (transition time, new state) pairs.
    The path sampling is conditional on the initial state
    but it is not conditional on the final state.
    So this is a 'forward' rather than a 'bridge' sampling.
    It does not require any matrix exponential computation.
    @param state_in: initial state
    @param blen_in: length of the branch
    @param rates: the rate away from each state
    @param P: transition matrix conditional on leaving a state
    """
    state = state_in
    nstates = len(rates)
    blen_accum = 0
    while True:
        rate = rates[state]
        scale = 1 / rate
        b = np.random.exponential(scale=scale)
        blen_accum += b
        if blen_accum >= blen_in:
            return
        distn = P[state]
        new_state = random_category(distn)
        yield blen_accum, new_state

def get_rejection_sample(
        total_length, initial_state, final_state, rates, P):
    """
    @param total_length: length of the path history in continuous time
    @param initial_state: state index at one end of the history
    @param final_state: state index at the other end of the history
    @param rates: rates away from the states
    @param P: substitution distributions conditional on instantaneous change
    """
    accepted_path = []
    while not accepted_path:
        length_state_pairs = list(gen_branch_history_sample(
            initial_state, total_length, rates, P))
        if length_state_pairs:
            obs_final_length, obs_final_state = length_state_pairs[-1]
            if obs_final_state == final_state:
                accum = 0
                state = initial_state
                for blen, next_state in length_state_pairs:
                    accepted_path.append((state, blen-accum))
                    state = next_state
                    accum = blen
                accepted_path.append((state, total_length-accum))
        elif initial_state == final_state:
            accepted_path.append((initial_state, total_length))
    return accepted_path


def build_single_history_table(conn, table, states, f_sample):
    """
    @param conn: database connection
    @param table: validated alphanumeric table name
    @param states: ordered list of integer states
    @param f_sample: rejection sampling function
    """

    # create the table
    cursor = conn.cursor()
    s = (
            'create table if not exists {table} ('
            'segment integer, '
            'state integer, '
            'blen real, '
            'primary key (segment))'
            ).format(table=table)
    cursor.execute(s)
    conn.commit()

    # populate the table
    accepted_path = f_sample()
    for segment_index, (state_index, blen) in enumerate(accepted_path):
        triple = (segment_index, states[state_index], blen)
        s = 'insert into {table} values (?, ?, ?)'.format(table=table)
        t = triple
        cursor.execute(s, t)
    conn.commit()


def main(args):

    # define the number of histories to sample
    nsamples = args.nsamples

    # validate table name or make a default table name
    if args.table is None:
        if nsamples == 1:
            table_name = 'history'
        else:
            table_name = 'histories'
    else:
        if not args.table.isalnum():
            raise Exception('table name must be alphanumeric')
        table_name = args.table

    # read the rate matrix
    conn = sqlite3.connect(args.rates)
    cursor = conn.cursor()
    cursor.execute('select source, sink, rate from rates')
    source_sink_rate = list(cursor)
    conn.close()

    # define the set of states and the rate matrix
    sources, sinks, rates = zip(*source_sink_rate)
    states = sorted(set(sources + sinks))
    nstates = len(states)
    s_to_i = dict((s, i) for i, s in enumerate(states))
    pre_Q = np.zeros((nstates, nstates), dtype=float)
    for si, sj, rate in source_sink_rate:
        pre_Q[s_to_i[si], s_to_i[sj]] = rate
    Q = pre_Q - np.diag(np.sum(pre_Q, axis=1))

    # do a conversion for rejection sampling
    rates, P = decompose_rates(Q)

    # partially evaluate the rejection sampling function
    f_sample = functools.partial(
            get_rejection_sample,
            args.elapsed,
            s_to_i[args.initial],
            s_to_i[args.final],
            rates, P)

    # sample some stuff
    conn = sqlite3.connect(args.outfile)
    build_single_history_table(conn, table_name, states, f_sample)
    conn.close()


if __name__ == '__main__':

    # define the ctmc bridge sampling methods
    method_names = (
            'rejection',
            #'modified-rejection',
            #'uniformization',
            #'direct',
            )

    # construct the command line
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--method', choices=method_names, default='rejection',
            help='use this sampling method')
    parser.add_argument('--initial', type=nonneg_int, required=True,
            help='initial state')
    parser.add_argument('--final', type=nonneg_int, required=True,
            help='final state')
    parser.add_argument('--elapsed', type=pos_float, default=1.0,
            help='final state')
    parser.add_argument('--rates', default='rate.matrix.db',
            help='input rate matrix as an sqlite3 database file')
    parser.add_argument('--nsamples', type=nonneg_int, default=4,
            help='sample this many endpoint-conditioned histories')
    parser.add_argument('--outfile', default='histories.db',
            help='output path samples as an sqlite3 database file')
    parser.add_argument('--table',
            help='name of table to create in the new database')
    main(parser.parse_args())

