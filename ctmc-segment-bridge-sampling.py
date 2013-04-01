"""
Sample the state history of a time interval given initial and final state.

For some endpoint conditioned path sampling algorithms, see
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2818752/
Also see the research report from the oxford summer school
for computational biology titled
conditional path sampling in continuous-time markov chains.
author last names are biswas, chan, xiong, tataru, herman, hobolth.
Actually the summer school pdf has errors in the equations,
so use the Hobolth and Stone 2009 paper instead.
"""

#XXX implement more path sampling methods
# * uniformization
# * ctmc uniformization

import functools
import argparse
import sqlite3
import math
import random

import numpy as np
import scipy.linalg
import scipy.optimize

import cmedbutil


#XXX this is copypasted
def gen_branch_history_sample(state_in, blen_in, rates, P, t0=0.0):
    """
    Forward path sampling along a branch with a known initial state.
    Yield (transition time, new state) pairs.
    The path sampling is conditional on the initial state
    but it is not conditional on the final state.
    So this is a 'forward' rather than a 'bridge' sampling.
    It does not require any matrix exponential computation.
    @param state_in: initial state
    @param blen_in: time until end of branch
    @param rates: the rate away from each state
    @param P: transition matrix conditional on leaving a state
    @param t0: initial time
    """
    t = t0
    state = state_in
    nstates = len(rates)
    while True:
        rate = rates[state]
        scale = 1 / rate
        delta_t = np.random.exponential(scale=scale)
        t += delta_t
        if t >= blen_in:
            return
        distn = P[state]
        state = cmedbutil.random_category(distn)
        yield t, state

def gen_modified_branch_history_sample(
        initial_state, final_state, blen_in, rates, P, t0=0.0):
    """
    This is a helper function for Nielsen modified rejection sampling.
    Yield (transition time, new state) pairs.
    The idea is to sample a path which may need to be rejected,
    and it is slightly clever in the sense that the path does not
    need to be rejected as often as do naive forward path samples.
    In more detail, this path sampler will generate paths
    conditional on at least one change occurring on the path,
    when appropriate.
    @param initial_state: initial state
    @param final_state: initial state
    @param blen_in: length of the branch
    @param rates: the rate away from each state
    @param P: transition matrix conditional on leaving a state
    @param t0: initial time
    """
    t = t0
    state = initial_state
    if state != final_state:
        rate = rates[initial_state]
        u = random.random()
        delta_t = -math.log1p(u*math.expm1(-blen_in*rate)) / rate
        t += delta_t
        if t >= blen_in:
            return
        distn = P[state]
        state = cmedbutil.random_category(distn)
        yield t, state
    for t, state in gen_branch_history_sample(state, blen_in, rates, P, t0=t):
        yield t, state


def get_naive_rejection_sample(
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
        t_state_pairs = list(gen_branch_history_sample(
            initial_state, total_length, rates, P))
        if t_state_pairs:
            obs_final_length, obs_final_state = t_state_pairs[-1]
            if obs_final_state == final_state:
                accum = 0
                state = initial_state
                for blen, next_state in t_state_pairs:
                    accepted_path.append((state, blen-accum))
                    state = next_state
                    accum = blen
                accepted_path.append((state, total_length-accum))
        elif initial_state == final_state:
            accepted_path.append((initial_state, total_length))
    return accepted_path

#XXX too much copypaste
def get_modified_rejection_sample(
        total_length, initial_state, final_state, rates, P):
    """
    If applicable, condition on at least one change.
    This modification often associated with Rasmus Nielsen (2002).
    @param total_length: length of the path history in continuous time
    @param initial_state: state index at one end of the history
    @param final_state: state index at the other end of the history
    @param rates: rates away from the states
    @param P: substitution distributions conditional on instantaneous change
    """
    accepted_path = []
    while not accepted_path:
        t_state_pairs = list(gen_modified_branch_history_sample(
            initial_state, final_state, total_length, rates, P))
        if t_state_pairs:
            obs_final_length, obs_final_state = t_state_pairs[-1]
            if obs_final_state == final_state:
                accum = 0
                state = initial_state
                for blen, next_state in t_state_pairs:
                    accepted_path.append((state, blen-accum))
                    state = next_state
                    accum = blen
                accepted_path.append((state, total_length-accum))
        elif initial_state == final_state:
            accepted_path.append((initial_state, total_length))
    return accepted_path

def sample_helper(cdf, u, t):
    """
    This is for sampling using an inverse transformation.
    It is mean to be used with functools.partial().
    @param cdf: returns a probability given t
    @param u: a target probability, e.g. drawn from the uniform unit interval
    @param t: a parameter value to be guessed
    @return: probability difference, for root finding
    """
    return cdf(t) - u

def conditional_waiting_time_cdf(
        initial_state, next_state, final_state,
        P, Q, w, U, U_inv, T, p_i, t):
    """
    This is a helper function for random sampling via root finding.
    It assumes that we have already sampled the next state,
    and that we need to sample the waiting time until that state.
    The endpoint conditioning makes this tricky.
    @param initial_state: initial state
    @param next_state: next state
    @param final_state: final state
    @param P: transition matrix depending on the time of the remaining branch
    @param Q: rate matrix
    @param w: eigenvalues
    @param U: part of eigendecomposition
    @param U_inv: part of eigendecomposition
    @param T: length of remaining branch
    @param p_i: a normalizing constant
    @param t: the waiting time
    """
    nstates = P.shape[0]
    qa = -Q[initial_state, initial_state]
    qai = Q[initial_state, next_state]
    Pab = P[initial_state, final_state]
    v = np.zeros(nstates, dtype=float)
    for j in range(nstates):
        denom = w[j] + qa
        if not denom:
            v[j] = t * math.exp(T * w[j])
        else:
            v[j] = -math.exp(T * w[j]) * math.expm1(-t*denom) / denom
    U_ij = U[next_state, :]
    U_inv_jb = U_inv[:, final_state]
    Fi = (qai / Pab) * np.sum(U_ij * U_inv_jb * v)
    return Fi / p_i


def gen_direct_sample(
        total_length, initial_state, final_state,
        Q, w, U, U_inv):
    """
    Use direct sampling as opposed to rejection sampling or uniformization.
    This uses an eigendecomposition of the rate matrix.
    The U and w are such that Q = U * diag(w) * U_inv in matrix notation.
    Note that this does not mean that U is assumed to be
    an orthogonal matrix.
    @param total_length: length of the path history in continuous time
    @param initial_state: state index at one end of the history
    @param final_state: state index at the other end of the history
    @param Q: rate matrix
    @param w: eigenvalues part of the rate matrix eigendecomposition
    @param U: square matrix part of the rate matrix eigendecomposition
    @param U_inv: square matrix part of the rate matrix eigendecomposition
    """
    t = 0
    T = total_length
    nstates = Q.shape[0]
    state = initial_state
    while True:
        P = scipy.linalg.expm(Q*T)
        q = -np.diag(Q)
        qa = q[state]
        pa = math.exp(-qa*T) / P[state, state]
        if state == final_state:
            if random.random() < pa:
                return
        J = np.zeros(nstates, dtype=float)
        for j in range(nstates):
            denom = w[j] + qa
            if not denom:
                J[j] = T * math.exp(T * w[j])
            else:
                J[j] = (math.exp(T * w[j]) - math.exp(-T*qa)) / denom
        distn = np.empty(nstates, dtype=float)
        for i in range(nstates):
            if i == state:
                distn[i] = pa
            else:
                q_ai = Q[state, i]
                U_ij = U[i, :]
                U_inv_jb = U_inv[:, final_state]
                J_aj = J
                denom = P[state, final_state]
                distn[i] = q_ai * np.sum(U_ij * U_inv_jb * J_aj) / denom
        """
        # check invariants
        x = sum(distn[j] for j in range(nstates) if j != state)
        if state == final_state:
            if not np.allclose(x, 1):
                raise Exception('distribution is wrong: %s %s' % (x, distn))
        else:
            if not np.allclose(x, 1 - pa):
                raise Exception('distribution is wrong')
        """
        #print distn
        #assert_stochastic_vector(distn)
        #raise Exception('it is OK')
        # sample the next state conditional on a state change
        change_distn = np.array(distn)
        change_distn[state] = 0
        change_distn /= np.sum(change_distn)
        next_state = cmedbutil.random_category(change_distn)
        #
        # use root-finding to get the waiting time
        #
        partial_cdf = functools.partial(
                conditional_waiting_time_cdf,
                state, next_state, final_state,
                P, Q, w, U, U_inv, T, distn[next_state])
        u = random.random()
        f = functools.partial(sample_helper, partial_cdf, u)
        delta_t = scipy.optimize.brentq(f, 0, T)
        t += delta_t
        if t > total_length:
            return
        T -= delta_t
        state = next_state
        yield t, state

def get_direct_sample(
        total_length, initial_state, final_state,
        Q, w, U, U_inv):
    accepted_path = []
    t_state_pairs = list(gen_direct_sample(
        total_length, initial_state, final_state, Q, w, U, U_inv))
    if t_state_pairs:
        obs_final_length, obs_final_state = t_state_pairs[-1]
        if obs_final_state == final_state:
            accum = 0
            state = initial_state
            for blen, next_state in t_state_pairs:
                accepted_path.append((state, blen-accum))
                state = next_state
                accum = blen
            accepted_path.append((state, total_length-accum))
    elif initial_state == final_state:
        accepted_path.append((initial_state, total_length))
    else:
        raise Exception('invalid path')
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
    rates, P = cmedbutil.decompose_rates(Q)

    # partially evaluate the rejection sampling function
    if args.method == 'naive-rejection':
        f_sample = functools.partial(
                get_naive_rejection_sample,
                args.elapsed,
                s_to_i[args.initial],
                s_to_i[args.final],
                rates, P)
    elif args.method == 'modified-rejection':
        f_sample = functools.partial(
                get_modified_rejection_sample,
                args.elapsed,
                s_to_i[args.initial],
                s_to_i[args.final],
                rates, P)
    elif args.method == 'direct':
        # FIXME this can be improved by symmetrizing
        w, U = scipy.linalg.eig(Q, right=True)
        U_inv = scipy.linalg.inv(U)
        if not np.allclose(np.eye(nstates), np.dot(U, U_inv)):
            raise Exception('bad eigenvectors')
        if not np.allclose(Q, np.dot(U, np.dot(np.diag(w), U_inv))):
            raise Exception('bad eigendecomposition')
        f_sample = functools.partial(
                get_direct_sample,
                args.elapsed,
                s_to_i[args.initial],
                s_to_i[args.final],
                Q, w, U, U_inv)
    else:
        raise Exception('the requested sampling method is not implemented')

    # sample some stuff
    conn = sqlite3.connect(args.outfile)
    if nsamples == 1:
        build_single_history_table(
                conn, table_name, states, f_sample)
    else:
        build_multiple_histories_table(
                conn, table_name, nsamples, states, f_sample)
    conn.close()


if __name__ == '__main__':

    # define the ctmc bridge sampling methods
    method_names = (
            'naive-rejection',
            'modified-rejection',
            'direct',
            #'uniformization',
            )

    # construct the command line
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--method', choices=method_names,
            default='naive-rejection',
            help='use this sampling method')
    parser.add_argument('--initial', type=cmedbutil.nonneg_int, required=True,
            help='initial state')
    parser.add_argument('--final', type=cmedbutil.nonneg_int, required=True,
            help='final state')
    parser.add_argument('--elapsed', type=cmedbutil.pos_float, default=1.0,
            help='elapsed time')
    parser.add_argument('--rates', default='rate.matrix.db',
            help='input rate matrix as an sqlite3 database file')
    parser.add_argument('--nsamples', type=cmedbutil.pos_int, default=4,
            help='sample this many endpoint-conditioned histories')
    parser.add_argument('--outfile', default='path.histories.db',
            help='output path samples as an sqlite3 database file')
    parser.add_argument('--table',
            help='name of table to create in the new database')
    main(parser.parse_args())

