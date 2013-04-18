"""
This script is purely for testing.

A path history is sampled according to a not particularly well defined model.
This script could possibly be generalized into an endpoint-conditioned
poisson process path sampling tool.
The output path sample format is just a single path sample.
The output tree sample format is a tree history sample
for a single history and for a single alignment column.
"""

import argparse
import random
import sqlite3

import numpy as np

import cmedbutil
from cmedbutil import pairwise


def write_path_history(sampled_segments, filename):
    """
    Write a file containing a sampled path history.
    @param sampled_segments: sequence of sampled (blen, state) pairs
    @param filename: write this sqlite3 database file
    """

    # open a database file for writing
    conn = sqlite3.connect(filename)
    cursor = conn.cursor()

    # initialize the output table
    cursor = conn.cursor()
    s = (
            'create table if not exists history ('
            'segment integer, '
            'state integer, '
            'blen real, '
            'primary key (segment))')
    cursor.execute(s)

    # populate the database
    for i, (blen, state) in enumerate(sampled_segments):
        segment = i
        s = 'insert into history values (?, ?, ?)'
        t = (segment, state, blen)
        cursor.execute(s, t)
    conn.commit()

    # close the output database
    conn.close()


def write_tree_history(sampled_segments, filename):
    """
    Write a file containing a sampled tree history.
    @param sampled_segments: sequence of sampled (blen, state) pairs
    @param filename: write this sqlite3 database file
    """

    # open a database file for writing
    conn = sqlite3.connect(filename)
    cursor = conn.cursor()

    # initialize the output table
    cursor = conn.cursor()
    s = (
            'create table if not exists histories ('
            'history integer, '
            'offset integer, '
            'segment integer, '
            'va integer, '
            'vb integer, '
            'blen real, '
            'state integer, '
            'primary key (history, offset, segment))')
    cursor.execute(s)

    # populate the database
    history_index = 0
    offset = 1
    for i, (blen, state) in enumerate(sampled_segments):
        segment = i
        va = i
        vb = i+1
        s = 'insert into histories values (?, ?, ?, ?, ?, ?, ?)'
        t = (history_index, offset, segment, va, vb, blen, state)
        cursor.execute(s, t)
    conn.commit()

    # close the output database
    conn.close()



def main(args):

    # check that we have at least three states
    if args.nstates < 3:
        raise Exception('not enough states')

    # check that we have at least three segments
    if args.nsegments < 3:
        raise Exception('not enough segments')
    
    # check that the initial and final states are in bounds
    if args.initial not in range(args.nstates):
        raise Exception('initial state is out of bounds')
    if args.final not in range(args.nstates):
        raise Exception('final state is out of bounds')

    # get the jump sequence
    jump_sequence = []
    full_state_set = set(range(args.nstates))
    for i in range(args.nsegments):
        if i == 0:
            state = args.initial
        elif i == args.nsegments - 1:
            state = args.final
        else:
            prev_state = jump_sequence[i-1]
            if i == args.nsegments - 2:
                taboo_set = {prev_state}
            else:
                taboo_set = {prev_state, args.final}
            state_choices = list(full_state_set - taboo_set)
            state = random.choice(state_choices)
        jump_sequence.append(state)

    # sample the nstates-1 substitution times along the path
    subs_times = np.random.uniform(0, args.blen, args.nsegments-1)

    # sample a few non-substitution times along the path
    non_subs_times = np.random.uniform(0, args.blen, 3)

    # get the sorted times, along with the binary substitution annotation
    annotated_times = sorted(
            [(t, 1) for t in subs_times] +
            [(t, 0) for t in non_subs_times] +
            [(0.0, 0), (args.blen, 0)])

    # get the segments with their lengths and their sampled states
    jump_index = 0
    sampled_segments = []
    for i, ((ta, sa), (tb, sb)) in enumerate(pairwise(annotated_times)):
        if sa:
            jump_index += 1
        segment = ((tb - ta), jump_sequence[jump_index])
        sampled_segments.append(segment)

    # construct a tree history file from the sampled path
    if args.tree_history:
        write_tree_history(sampled_segments, args.tree_history)

    # construct a path history file from the sampled path
    if args.path_history:
        write_path_history(sampled_segments, args.path_history)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)

    # initial state, final state, elapsed time
    parser.add_argument('--initial', default=0, type=cmedbutil.nonneg_int,
            help='initial state')
    parser.add_argument('--final', default=1, type=cmedbutil.nonneg_int,
            help='final state')
    parser.add_argument('--blen', default=1, type=cmedbutil.pos_float,
            help='branch length associated with the path')

    # other cmdline arguments defining the process and conditions
    parser.add_argument('--nstates', default=4, type=cmedbutil.pos_int,
            help='total number of states in the process')
    parser.add_argument('--nsegments', default=6, type=cmedbutil.pos_int,
            help='number of segments')

    # output paths
    parser.add_argument('--path-history', default='toy.path.history.db',
            help='output path history as a sqlite3 file')
    parser.add_argument('--tree-history', default='toy.tree.history.db',
            help='output tree history as a sqlite3 file')
    main(parser.parse_args())

