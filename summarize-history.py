"""
Summarize sampled histories by reporting averages of summary statistics.

The output of this script should be an sqlite3 database
that has the same format as the output of the script
that computes the exact endpoint-conditioned expectations.
The difference is that this script does not average
over all possible histories.
Instead it averages over only a small set of sampled histories.
.
Maybe later allow only a single history.
.
Input schema.
[multiple path (not tree) endpoint-conditioned history samples]
table histories (possibly you can accept this as input)
history integer, (this is just a part of identification of the entry)
segment integer, (this is another part of identification of the entry)
state integer, (this is part of the guts of the path)
blen real (this is another part of the guts of the path)
.
Output schema.
table wait
initial integer, final integer, state integer, wait real
table usage
initial integer, final integer, source integer, sink integer, usage real
"""

import argparse
import sqlite3
import math
import itertools

import numpy as np
import scipy.linalg


#XXX this could go into a module of recipes
def pairwise(iterable):
    """
    This is an itertools recipe.
    s -> (s0,s1), (s1,s2), (s2, s3), ...
    """
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)

#XXX this should go into a separate module
def nonneg_int(x):
    x = int(x)
    if x < 0:
        raise argparse.ArgumentTypeError(
                'value must be a non-negative integer')
    return x

#XXX this should go into a separate module
def pos_int(x):
    x = int(x)
    if x < 1:
        raise argparse.ArgumentTypeError(
                'value must be a positive integer')
    return x

#XXX this should go into a separate module
def pos_float(x):
    x = float(x)
    if x <= 0:
        raise argparse.ArgumentTypeError(
                'value must be a positive floating point number')
    return x

#XXX this should go into a separate module
def assert_stochastic_vector(v):
    if np.any(v < 0) or np.any(1 < v):
        raise Exception(
                'entries of a finite distribution vector should be in '
                'the inclusive interval [0, 1]')
    if not np.allclose(np.sum(v), 1):
        raise Exception(
                'entries of a finite distribution vector should sum to 1')

#XXX this should go into a separate module
def assert_rate_matrix(Q):
    if not np.allclose(np.sum(Q, axis=1), 0):
        raise Exception('expected rate matrix rows to sum to zero')
    if np.any(np.diag(Q) > 0):
        raise Exception('expected rate matrix diagonals to be non-positive')
    if np.any(Q - np.diag(np.diag(Q)) < 0):
        raise Exception('expected rate matrix off-diagonals to be non-negative')

#XXX this should go into a separate module
def assert_equilibrium(Q, distn):
    if not np.allclose(np.dot(distn, Q), 0):
        raise Exception('the distribution is not at equilibrium')

#XXX this should go into a separate module
def assert_detailed_balance(Q, distn):
    S = (Q.T * distn).T
    if not np.allclose(S, S.T):
        raise Exception('the detailed balance equations are not met')


def main(args):

    # read and validate the rate matrix info from the sqlite3 database file
    conn = sqlite3.connect(args.rates)
    cursor = conn.cursor()
    states, distn, Q = get_rate_matrix_info(cursor)
    conn.close()

    # compute the map from state to state index
    s_to_i = dict((s, i) for i, s in enumerate(states))

    # define the initial states of interest
    if args.initial is None:
        initial_states = states
    else:
        initial_state = args.initial
        if initial_state not in states:
            raise Exception('unknown initial state %s' % initial_state)
        initial_states = [initial_state]

    # define the final states of interest
    if args.final is None:
        final_states = states
    else:
        final_state = args.final
        if final_state not in states:
            raise Exception('unknown final state %s' % final_state)
        final_states = [final_state]

    # extract the amount of time along the path
    T = args.elapsed

    # Compute the matrix exponential M.
    # This matrix is often represented as P but the Holmes-Rubin (2002)
    # paper uses the notation M.
    M = scipy.linalg.expm(Q*T)

    # Use the properties of time-reversibility to symmetrize the rate matrix.
    # The symmetric matrix S uses the Holmes-Rubin (2002) notation.
    r = np.sqrt(distn)
    S = (Q.T * r).T / r
    if not np.allclose(S, S.T):
        raise Exception('internal error constructing symmetrized matrix')

    # compute the symmetric eigendecomposition of the symmetric matrix
    w, V = scipy.linalg.eigh(S)

    # compute the exact expectations of the summary statistics
    J = get_decay_mode_interactions(w, T)

    # create the output database file and initialize the cursor
    conn = sqlite3.connect(args.outfile)
    cursor = conn.cursor()

    # define the table of waiting times
    s = (
            'create table if not exists wait ('
            'initial integer, '
            'final integer, '
            'state integer, '
            'wait real, '
            'primary key (initial, final, state))')
    cursor.execute(s)
    conn.commit()

    # define the table of transition usages
    s = (
            'create table if not exists usage ('
            'initial integer, '
            'final integer, '
            'source integer, '
            'sink integer, '
            'usage real, '
            'primary key (initial, final, source, sink))')
    cursor.execute(s)
    conn.commit()

    # populate the waiting times table
    for initial_state in initial_states:
        for final_state in final_states:
            a = s_to_i[initial_state]
            b = s_to_i[final_state]
            Tau = get_intermediate_matrix(a, b, distn, V, J)
            wait_times = get_expected_wait_time(a, b, M, Tau)
            if not np.allclose(np.sum(wait_times), T):
                raise Exception(
                        'waiting time expectations '
                        'should add up to the elapsed time')
            if np.any(wait_times < 0):
                raise Exception(
                        'waiting time expectations '
                        'should be non-negative')
            for i, wait in enumerate(wait_times):
                state = s_to_i[i]
                s = 'insert into wait values (?, ?, ?, ?)'
                t = (initial_state, final_state, state, wait)
                cursor.execute(s, t)
    conn.commit()

    # populate the transition usage expectation table
    for initial_state in initial_states:
        for final_state in final_states:
            a = s_to_i[initial_state]
            b = s_to_i[final_state]
            Tau = get_intermediate_matrix(a, b, distn, V, J)
            usages = get_expected_transition_usage(a, b, Q, M, Tau)
            if np.any(usages < 0):
                raise Exception(
                        'transition usage count expectations '
                        'should be non-negative')
            for i, source_state in enumerate(states):
                for j, sink_state in enumerate(states):
                    if i != j:
                        usage = usages[i, j]
                        s = 'insert into usage values (?, ?, ?, ?, ?)'
                        t = (
                                initial_state, final_state,
                                source_state, sink_state,
                                usage)
                        cursor.execute(s, t)
    conn.commit()

    # close the output database connection
    conn.close()


def get_summary_from_history(states, history):
    """
    Each history is a sample path.
    Each sample path is a sequence of (state, wait) pairs.
    The output wait times is a vector of expected wait times per history.
    The output transition time matrix is a square ndarray
    of expected transition counts per history.
    Because this function acts on only a single history,
    the expectations are just the counts of observed values.
    @param states: ordered states
    @param history: a single sample path
    @return: wait_times, transition_counts
    """

    # Precompute the map from states to state indices.
    s_to_i = dict((s, i) for i, s in enumerate(states))
    nstates = len(states)

    # Get the wait time expectations from the sampled paths.
    wait_times = np.zeros(nstates, dtype=float)
    for state, wait in history:
        wait_times[state] += wait

    # Count the number of each transition along the path.
    transition_counts = np.zeros((nstates, nstates), dtype=float)
    for ((state_a, wait_a), (state_b, wait_b)) in pairwise(history):
        a = s_to_i[state_a]
        b = s_to_i[state_b]
        transition_counts[a, b] += 1

    # Return the wait times and transition counts for this single history.
    return wait_times, transition_counts


def get_summary_expectation_from_histories(states, histories):
    """
    Each history is a sample path.
    Each sample path is a sequence of (state, blen) pairs.
    The output wait times is a vector of expected wait times per history.
    The output transition time matrix is a square ndarray
    of expected transition counts per history.
    @param states: ordered states
    @param histories: sequence of sample paths
    @return: wait times expectation, transition count expectation
    """
    nhistories = len(histories)
    wait_time_expectation = np.zeros(nstates, dtype=float)
    transition_count_expectation = np.zeros((nstates, nstates), dtype=float)
    for history in histories:
        waits, trans = get_expectation_from_history(states, history)
        wait_time_expectation += waits
        transition_count_expectation += trans
    wait_time_expectation /= float(nhistories)
    transition_count_expectation /= float(nhistories)
    return wait_time_expectation, transition_count_expectation


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outfile', default='averages.db',
            help='output averages as an sqlite3 database file')
    main(parser.parse_args())

