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


def main(args):

    # read the path histories from the sqlite3 database file
    conn = sqlite3.connect(args.histories)
    cursor = conn.cursor()
    cursor.execute('select history, segment, state, blen from histories')
    data = sorted(cursor)
    states = sorted(cursor.execute('select state from histories'))
    conn.close()

    # Define the map from state to rate matrix index,
    # and count the number of different states in the history.
    s_to_i = dict((s, i) for i, s in enumerate(states))
    nstates = len(states)

    # parse the sample path histories from the table
    histories = []
    segment = []
    for history_index, segment_index, state, blen in data:
        if history_index == len(histories):
            histories.append(segment)
            segment = []
        if history_index != len(histories) - 1:
            raise Exception('invalid history index')
        if segment_index != len(segment):
            raise Exception('invalid segment index')
        segment.append((state, blen))
    histories.append(segment)

    # get the averages from the histories
    average_weight_times, average_transition_counts = get_summary_averages(
            states, histories)

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


#FIXME this needs the endpoint conditioning!
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
        i = s_to_i[state_i]
        j = s_to_i[state_j]
        transition_counts[a, b] += 1

    # Return the wait times and transition counts for this single history.
    return wait_times, transition_counts


#FIXME this needs the endpoint conditioning!
def get_summary_averages(states, histories):
    """
    Each history is a sample path.
    Each sample path is a sequence of (state, blen) pairs.
    The output wait times is a vector of expected wait times per history.
    The output transition time matrix is a square ndarray
    of expected transition counts per history.
    @param states: ordered states
    @param histories: sequence of sample paths
    @return: wait times averages, transition count averages
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
    parser.add_argument('--histories', default='histories.db',
            help='input path histories as an sqlite3 database file')
    parser.add_argument('--outfile', default='averages.db',
            help='output averages as an sqlite3 database file')
    main(parser.parse_args())

