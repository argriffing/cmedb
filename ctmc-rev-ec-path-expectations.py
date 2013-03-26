"""
Compute endpoint-conditioned expectations.

The meanings of the components of the script name are as follows.
This is related to continuous-time Markov chain processes hence ctmc.
The processes are assumed to be statistically time-reversible hence rev.
The ec means endpoint-conditioned,
which is sometimes also referred to as bridged, like a Brownian bridge
or bridge sampling.
The process is assumed to act in continuous time along a path rather
than along a star tree or along a general tree,
hence the word 'path' in the script name.
Finally, we want to compute conditional expectations of some stuff.
The conditional expectations are the quantities in equation (4)
of Holmes and Rubin (2002).
These are the expected number of paths that start in state i,
the expected wait in state i (i.e. the amount of time spent in i)
and the expected usage of transition i->j.
.
The bigger picture is that this script is about computing
expected summary statistics of path histories given the endpoints
and the process and the path length for the purpose of
testing path sampling algorithms by comparing averages
of their summary statistics vs. the theoretical
expected values of the summary statistics.
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

def get_decay_mode_interactions(w, T):
    """
    This returns the J ndarray.
    @param w: eigenvalues
    @param T: elapsed time
    @return: a matrix of interactions between pairs of decay modes
    """
    nstates = w.shape[0]
    J = np.empty((nstates, nstates))
    v = np.exp(w * T)
    for k in range(nstates):
        for l in range(nstates):
            denom = w[k] - w[l]
            if denom:
                J[k, l] = (v[k] - v[l]) / denom
            else:
                J[k, l] = T * v[k]
    return J

def get_intermediate_matrix(a, b, distn, V, J):
    """
    This returns the Tau ndarray used by the get_expected_* functions.
    @param a: initial endpoint state
    @param b: final endpoint state
    @param distn: prior equilibrium distribution of the time-reversible process
    @param V: orthogonal matrix as an ndarray
    @param J: decay mode interactions including eigenvalue information
    @return: an intermediate matrix for computing expectations
    """
    #FIXME: this is translated from (4) in Holmes-Rubin (2002)
    # and it could be rewritten in linear algebra notation.
    # presumably Asger Hobolth has done this using R.
    nstates = V.shape[0]
    tau_prefix = np.empty_like(J)
    for i in range(nstates):
        for j in range(nstates):
            num = distn[i] * distn[b]
            den = distn[a] * distn[j]
            tau_prefix[i, j] = np.sqrt(num / den)
    tau_suffix = np.zeros_like(J)
    for i in range(nstates):
        for j in range(nstates):
            for k in range(nstates):
                inner = np.sum(V[j, :] * V[b, :] * J[k, :])
                tau_suffix[i, j] += V[a, k] * V[i, k] * inner
    Tau = tau_prefix * tau_suffix
    return Tau

def get_expected_wait_time(a, b, M, Tau):
    """
    This is w hat in Eq. (4) of Holmes-Rubin 2002.
    @param a: initial endpoint state
    @param b: final endpoint state
    @param M: matrix exponential expm(T*R)
    @param Tau: an intermediate matrix for computing expectations
    @return: a vector of expected wait times (amount of time spent in state i)
    """
    return np.diag(Tau) / M[a, b]

def get_expected_transition_usage(a, b, R, M, Tau):
    """
    This is u hat in Eq. (4) of Holmes-Rubin 2002.
    The diagonal entries of the returned ndarray are meaningless,
    so they are arbitrarily filled with zeros.
    @param a: initial endpoint state
    @param b: final endpoint state
    @param R: instantaneous transition rate matrix
    @param M: matrix exponential expm(T*R)
    @param Tau: an intermediate matrix for computing expectations
    @return: ndarray with expected count of each transition
    """
    U = (R * Tau) / M[a, b]
    np.fill_diagonal(U, 0)
    return U

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


def build_multiple_histories_table(conn, table, nsamples, states, f_sample):
    """
    @param conn: database connection
    @param table: validated alphanumeric table name
    @param nsamples: sample this many path histories
    @param states: ordered list of integer states
    @param f_sample: rejection sampling function
    """

    # create the table
    cursor = conn.cursor()
    s = (
            'create table if not exists {table} ('
            'history integer, '
            'segment integer, '
            'state integer, '
            'blen real, '
            'primary key (history, segment))'
            ).format(table=table)
    cursor.execute(s)
    conn.commit()

    # populate the table
    for history_index in range(nsamples):
        accepted_path = f_sample()
        for segment_index, (state_index, blen) in enumerate(accepted_path):
            s = 'insert into {table} values (?, ?, ?, ?)'.format(table=table)
            t = (history_index, segment_index, states[state_index], blen)
            cursor.execute(s, t)
    conn.commit()


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
    assert_stochastic_vector(distn)

    # assert that the rate matrix is actually a rate matrix
    assert_rate_matrix(Q)

    # assert that the distribution is at equilibrium w.r.t. the rate matrix
    assert_equilibrium(Q, distn)

    # assert that the detailed balance equations are met
    assert_detailed_balance(Q, distn)

    # return the validated inputs describing the stochastic process
    return states, distn, Q


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

    # Input requirements for Holmes-Rubin (2002) expectations:
    # - time-reversible rate matrix
    # - equilibrium distribution
    # - branch length
    #
    # Should the initial and final states be specified,
    # or should they be left free?
    # Possibly allow them to be optionally specified,
    # and compute all combinations of unspecified endpoints.
    #
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--initial', type=nonneg_int,
            help='initial state')
    parser.add_argument('--final', type=nonneg_int,
            help='final state')
    parser.add_argument('--elapsed', type=pos_float, default=1.0,
            help='elapsed time')
    parser.add_argument('--rates', default='rate.matrix.db',
            help='time-reversible rate matrix as an sqlite3 database file')
    parser.add_argument('--outfile', default='expectations.db',
            help='output expectations as an sqlite3 database file')
    main(parser.parse_args())

