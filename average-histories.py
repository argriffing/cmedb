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


#XXX this could go into a module of recipes
def pairwise(iterable):
    """
    This is an itertools recipe.
    s -> (s0,s1), (s1,s2), (s2, s3), ...
    """
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)


def main(args):

    # read the path histories from the sqlite3 database file
    conn = sqlite3.connect(args.histories)
    cursor = conn.cursor()
    cursor.execute('select history, segment, state, blen from histories')
    data = sorted(cursor)
    cursor.execute('select state from histories')
    states = sorted(set(t[0] for t in cursor))
    conn.close()

    # Define the map from state to rate matrix index,
    # and count the number of different states in the history.
    s_to_i = dict((s, i) for i, s in enumerate(states))
    nstates = len(states)

    # parse the sample path histories from the table
    histories = []
    segment = []
    for row in data:
        history_index, segment_index, state, blen = row
        if history_index > len(histories):
            histories.append(segment)
            segment = []
        if history_index != len(histories):
            print row
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
    for initial_state in states:
        for final_state in states:
            a = s_to_i[initial_state]
            b = s_to_i[final_state]
            wait_times = average_weight_times[a, b]
            if np.any(wait_times < 0):
                raise Exception(
                        'waiting time expectations '
                        'should be non-negative')
            for state, wait in zip(states, wait_times):
                s = 'insert into wait values (?, ?, ?, ?)'
                t = (initial_state, final_state, state, wait)
                cursor.execute(s, t)
    conn.commit()

    # populate the transition usage expectation table
    for initial_state in states:
        for final_state in states:
            a = s_to_i[initial_state]
            b = s_to_i[final_state]
            usages = average_transition_counts[a, b]
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
    n = len(states)

    # Get the wait time expectations from the sampled paths.
    wait_times = np.zeros(n, dtype=float)
    for state, wait in history:
        i = s_to_i[state_a]
        wait_times[i] += wait

    # Count the number of each transition along the path.
    transition_counts = np.zeros((n, n), dtype=float)
    for ((state_a, wait_a), (state_b, wait_b)) in pairwise(history):
        i = s_to_i[state_a]
        j = s_to_i[state_b]
        transition_counts[i, j] += 1

    # Return the wait times and transition counts for this single history.
    return wait_times, transition_counts


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

    # Precompute the map from states to state indices.
    s_to_i = dict((s, i) for i, s in enumerate(states))
    n = len(states)

    # Compute the summary statistics of the sampled paths.
    nhistories = len(histories)
    average_wait_times = np.zeros((n, n, n), dtype=float)
    average_transition_counts = np.zeros((n, n, n, n), dtype=float)
    for history in histories:
        initial_state, initial_blen = history[0]
        final_state, final_blen = history[-1]
        a = s_to_i[initial_state]
        b = s_to_i[final_state]
        waits, trans = get_summary_from_history(states, history)
        average_wait_times[a, b] += waits
        average_transition_counts[a, b] += trans
    average_wait_times /= float(nhistories)
    average_transition_counts /= float(nhistories)
    return average_wait_times, average_transition_counts


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--histories', default='histories.db',
            help='input path histories as an sqlite3 database file')
    parser.add_argument('--outfile', default='averages.db',
            help='output averages as an sqlite3 database file')
    main(parser.parse_args())

