"""
Sample a complicated Markov-modulated process along a path.

It could take a rate matrix and a state partition
and a blink birth rate and death rate and an elapsed time as inputs.
The output would be a couple of tables in an sqlite3 database.
The history table would map path segment indices
to (duration, primary state) pairs.
The blink table would map (part, path segment index) pairs
to (duration, blink state) pairs.
.
This script should not need to use expm.
Sampling on trees is not yet supported,
and sampling of multiple independent path histories is not yet supported.
Another future feature could be to specify the initial state
instead of drawing randomly from the stationary distribution.
Also, it could become useful to specify that some information about the
sampled history would be left out;
Perhaps all of the blink state would be known initially,
and a single blink state is known at the other endpoint
(it is known that the blink state corresponding to the
part of the final primary state is in the 'on' blink state),
but no information about the blink states inside of the path
is directly observable.
.
Possibly use networkx directed weighted graphs
as a representation of the sparse rate matrix.
"""

#FIXME redo this to not use numpy,
# and to more flexibly support sparse rate matrices
# without worrying about extra conversions between states
# and indices into ndarrays.
# This will mean not using the usual rate matrix reader.

# NOTE: numpy is used only for allclose() comparison

import argparse
import sqlite3
from collections import defaultdict
from itertools import permutations

import numpy as np
import networkx as nx

import cmedbutil

def get_moves(compound_state, dg, partition, rate_on, rate_off):
    """
    Get the successor states and their corresponding rates.
    Compound states have two substates.
    The first substate is the primary state which is an integer.
    The second substate is a map from partition part to a binary integer.
    @param compound_state: current compound state
    @param dg: directed graph of primary state transition rates
    @param partition: map from primary state to partition
    @param rate_on: a blinking rate
    @param rate_off: a blinking rate
    @return: a sequence of (successor compound state, rate) pairs
    """

    # get the set of indices of parts
    parts = set(partition.values())

    # expand the current compound state into its substates
    primary_state, blink_state = compound_state
    curr_part = partition[primary_state]

    # assert that the current part is blinked on
    if not blink_state[curr_part]:
        raise Exception('invalid compound state')

    # Compute the allowed moves out of the current compound state.
    # This is the sum of allowed primary transitions
    # and allowed blink toggles.
    move_rate_pairs = []

    # Add the allowed primary transitions.
    for b in dg.successors(primary_state):
        next_part = partition[b]
        if blink_state[next_part]:
            move = (b, dict(blink_state))
            rate = dg[primary_state][b]['weight']
            move_rate_pairs.append((move, rate))

    # Add the allowed blink rates.
    for blink_part in parts:
        if blink_part != curr_part:
            next_blink_state = dict(blink_state)
            if blink_state[blink_part] == 1:
                next_blink_state[blink_part] = 0
                rate = rate_off
            elif blink_state[blink_part] == 0:
                next_blink_state[blink_part] = 1
                rate = rate_on
            else:
                raise Exception('blink state of each part must be binary')
            move = (primary_state, next_blink_state)
            move_rate_pairs.append((move, rate))

    # Return the allowed moves out of the current state,
    # and their corresponding rates.
    return move_rate_pairs


def gen_branch_history_sample(
        primary_state_in, blink_state_in, blen_in,
        dg,
        partition, rate_on, rate_off,
        ):
    """
    This function is defined in analogy to gen_branch_history_sample.
    Path sampling along a branch with a known initial state.
    Yield (transition time, new primary state, new blink state) tuples.
    This function takes arguments divided into three groups.
    The first group defines the initial state and the path duration.
    The second group defines the primary process without regard to blinking;
    this has a weighted directed networkx graph for sparse rates.
    The third group defines the partition of the primary state space
    and the toggling rates of the modulating blinking processes.
    @param primary_state_in: initial integer state of the observable process
    @param blink_state_in: initial binary tuple state of the blinking process
    @param blen_in: length of the branch
    @param dg: rate matrix without regard to blinking
    @param partition: maps a primary state to a partition index for blinking
    @param rate_on: instantaneous off-to-on blinking rate
    @param rate_off: instantaneous on-to-off blinking rate
    """

    # The state space of this process is somewhat complicated.
    # The primary state is an integer.
    # The blink state is a map from partition part to binary blink state.
    # An invariant of the compound state
    # is that the blink state of part k must be 1 if the partition part of
    # the primary state is k.
    # Otherwise, any combination of primary states and blink states is allowed.
    primary_state = primary_state_in
    blink_state = dict(blink_state_in)
    blen_accum = 0
    while True:

        # Get the list of the allowed moves out of the current state,
        # and the corresponding rates associated with these moves.
        compound_state = (primary_state, blink_state)
        moves = get_moves(compound_state, dg, partition, rate_on, rate_off)
        successors, rates = zip(*moves)

        # Compute the total rate out of the current compound state,
        # and draw a random wait time that depends on this total rate.
        # If this wait time puts us over the allotted time period,
        # then we are done sampling the histories.
        total_rate = sum(rates)
        scale = 1.0 / total_rate
        blen_delta = np.random.exponential(scale=scale)
        blen_accum += blen_delta
        if blen_accum >= blen_in:
            return

        # Next we randomly pick a successor state
        # according to the rate proportions.
        pre_distn = np.array(rates, dtype=float)
        distn = pre_distn / np.sum(pre_distn)
        next_compound_state = successors[cmedbutil.random_category(distn)]
        primary_state, blink_state = next_compound_state

        # Yield the cumulative time and the new primary and blink states
        yield blen_accum, primary_state, blink_state


def get_sparse_rate_matrix_info(cursor):
    """
    This is a non-numpy customization of the usual get_rate_matrix_info.
    In this function we ignore the 'states' table with the state names,
    but we care about the sparse rates and the equilibrium distribution.
    Return a couple of things.
    The first thing is a map from a state to an equilibrium probability.
    The second thing is a sparse rate matrix as a networkx weighted digraph.
    @param cursor: sqlite3 database cursor
    @return: eq distn, rate graph
    """

    # get a set of all states
    cursor.execute(
            'select state from distn '
            'union '
            'select source from rates '
            'union '
            'select sink from rates '
            'union '
            'select state from states '
            )
    states = set(t[0] for t in cursor)
    nstates = len(states)

    # get the sparse equilibrium distribution
    cursor.execute('select state, prob from distn')
    distn = dict(cursor)

    # construct the rate matrix as a networkx weighted directed graph
    dg = nx.DiGraph()
    cursor.execute('select source, sink, rate from rates')
    for a, b, weight in cursor:
        dg.add_edge(a, b, weight=weight)

    # assert that the distribution has the right form
    if not all(0 <= p <= 1 for p in distn.values()):
        raise Exception(
                'equilibrium probabilities '
                'should be in the interval [0, 1]')
    if not np.allclose(sum(distn.values()), 1):
        raise Exception('equilibrium probabilities should sum to 1')

    # assert that rates are not negative
    if any(data['weight'] < 0 for a, b, data in dg.edges(data=True)):
        raise Exception('rates should be non-negative')

    # assert detailed balance
    for a, b in permutations(states, 2):
        if not np.allclose(
                distn[a] * dg[a][b]['weight'],
                distn[b] * dg[b][a]['weight'],
                ):
            raise Exception('expected detailed balance')

    # return the eq distn and the rate graph
    return distn, dg



def main(args):

    # read the sparse rate matrix from a database file
    conn = sqlite3.connect(args.rates)
    cursor = conn.cursor()
    distn, dg = get_sparse_rate_matrix_info(cursor)
    conn.close()

    # read the partition from a database file
    conn = sqlite3.connect(args.partition)
    cursor = conn.cursor()
    partition = dict(cursor.execute('select state, part from partition'))
    conn.close()

    # get the set of indices of parts
    parts = set(partition.values())
    nparts = len(parts)

    # Sample the initial state.
    # Do this by sampling the primary state according to its distribution,
    # then setting the corresponding blinking state to 'on',
    # then sampling the remaining blinking states
    # according to the bernoulli probability defined by the
    # blinking-on and blinking-off rates.
    states, probs = zip(*distn.items())
    primary_state_in = states[cmedbutil.random_category(probs)]
    p_blink_is_on = args.rate_on / (args.rate_on + args.rate_off)
    iid_blinks = [int(x) for x in np.random.binomial(1, p_blink_is_on, nparts)]
    blink_state_in = dict(zip(parts, iid_blinks))
    blink_state_in[partition[primary_state_in]] = 1

    # For the primary state and for each blink substate,
    # record the previous transition time or 0,
    # and also record its current state,
    # and also record its segment index.
    primary_state = primary_state_in
    primary_state_tm = 0.0
    primary_state_seg = 0
    blink_state = dict(blink_state_in)
    blink_state_tm = dict((part, 0.0) for part in parts)
    blink_state_seg = dict((part, 0) for part in parts)

    # sample the history
    primary_seg_state_duration_list = []
    blink_part_seg_state_duration_list = []
    for blen_accum, next_primary, next_blink in gen_branch_history_sample(
            primary_state_in, blink_state_in, args.duration,
            dg,
            partition, args.rate_on, args.rate_off,
            ):
        # Each time a transition happens,
        # this is an endpoint of one segment,
        # so add the segment to the appropriate history.
        if next_primary != primary_state:
            primary_seg_state_duration_list.append((
                primary_state_seg,
                primary_state,
                blen_accum - primary_state_tm))
            primary_state = next_primary
            primary_state_tm = blen_accum
            primary_state_seg += 1
        else:
            for part in parts:
                if next_blink[part] != blink_state[part]:
                    blink_part_seg_state_duration_list.append((
                        part,
                        blink_state_seg[part],
                        blink_state[part],
                        blen_accum - blink_state_tm[part]))
                    blink_state[part] = next_blink[part]
                    blink_state_tm[part] = blen_accum
                    blink_state_seg[part] += 1
    
    # Add trailing segments.
    primary_seg_state_duration_list.append((
        primary_state_seg,
        primary_state,
        args.duration - primary_state_tm))
    for part in parts:
        blink_part_seg_state_duration_list.append((
            part,
            blink_state_seg[part],
            blink_state[part],
            args.duration - blink_state_tm[part]))

    # open the output history database for writing
    conn = sqlite3.connect(args.outfile)
    cursor = conn.cursor()

    # initialize the primary state history table
    s = (
            'create table if not exists primary_history ('
            'segment integer, '
            'state integer, '
            'duration real, '
            'primary key (segment))')
    cursor.execute(s)
    conn.commit()

    # initialize the blink state history table
    s = (
            'create table if not exists blink_history ('
            'part integer, '
            'segment integer, '
            'state integer, '
            'duration real, '
            'primary key (part, segment))')
    cursor.execute(s)
    conn.commit()

    # populate the primary state history table
    for t in primary_seg_state_duration_list:
        s = 'insert into primary_history values (?, ?, ?)'
        cursor.execute(s, t)
    conn.commit()

    # populate the blink state history table
    for t in blink_part_seg_state_duration_list:
        s = 'insert into blink_history values (?, ?, ?, ?)'
        cursor.execute(s, t)
    conn.commit()

    # close the database file
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rate-on',
            type=cmedbutil.pos_float, default=1.0,
            help='rate at which blink states change from off to on')
    parser.add_argument('--rate-off',
            type=cmedbutil.pos_float, default=0.5,
            help='rate at which blink states change from on to off')
    parser.add_argument('--duration',
            type=cmedbutil.pos_float, default=1.0,
            help='total elapsed time along the path')
    parser.add_argument('--rates', default='rate.matrix.db',
            help=('continuous-time Markov substitution process '
                'as an sqlite3 database file'))
    parser.add_argument('--partition', default='partition.db',
            help=('a partition of the primary states '
                'as an sqlite3 database file'))
    parser.add_argument('-o', '--outfile', default='blink.path.db',
            help='create this sqlite3 database file')
    main(parser.parse_args())

