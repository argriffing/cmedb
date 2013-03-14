"""
Sample state histories along a tree.

This is unfinished.
The number of sampled histories is a command line option.
The root (for the purposes of sampling) is also given on the command line.
If the rate matrix is time-reversible then the choice of root does not matter.
The tree topology and branch lengths are given by a database.
The rates and initial state distribution are given by another database.
history integer
segment integer
va integer
vb integer
blen real
state integer
primary key (history, segment)
"""

import argparse
import sqlite3

import numpy as np


def main(args):

    # define the number of histories to sample
    nsamples = args.nsamples

    # open the rate matrix db file
    # read the rates
    # read the state distribution at the root
    conn = sqlite3.connect(args.rates)
    cursor = conn.cursor()
    cursor.execute(
            'select state from distn '
            'union '
            'select source from rates '
            'union '
            'select sink from rates '
            'union '
            'select state from states '
            )
    states = [t[0] for t in cursor]
    nstates = len(states)
    # define the map from state to rate matrix index
    s_to_i = dict((s, i) for i, s in enumerate(states))
    # construct the rate matrix
    cursor.execute('select source, sink, rate from rates')
    pre_Q = np.zeros((nstates, nstates), dtype=float)
    for si, sj, rate in cursor:
        pre_Q[s_to_i[si], s_to_i[sj]] = rate
    Q = pre_Q - np.diag(np.sum(pre_Q, axis=1))
    conn.close()

    # open the tree db file
    # read the association of node index pairs to edge indices
    # read the association of edge indices to branch lengths
    conn = sqlite3.connect(args.tree)
    cursor = conn.cursor()
    cursor.execute(
            'select va from topo '
            'union '
            'select vb from topo '
            )
    vertices = list(cursor)
    nvertices = len(vertices)
    cursor.execute('select edge, va, vb from topo')
    edge_va_vb_list = list(cursor)
    cursor.execute('select edge, blen from blen')
    edge_to_blen = dict(cursor)
    conn.close()

    # define the root node index for sampling;
    # the root is either user supplied,
    # or is taken to be the smallest node index otherwise
    if args.root is None:
        root = vertices[0]
    elif args.root in vertices:
        root = args.root
    else:
        raise Exception('the specified root should be a node index in the tree')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tree', required=True,
            help='unrooted tree as an sqlite3 database file')
    parser.add_argument('--rates', required=True,
            help=('continuous-time Markov substitution process '
                'as an sqlite3 database file'))
    parser.add_argument('--root', type=int,
            help='root node index')
    parser.add_argument('--nsamples', type=int, default=10,
            help='sample this many histories')
    main(parser.parse_args())

