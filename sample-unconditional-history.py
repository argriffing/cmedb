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
import networkx as nx


def sample_history(root, G_dag_in, edge_to_blen_in, distn, Q):
    """
    @param root: root of the tree
    @param G_dag_in: tree as a directed acyclic graph
    @param edge_to_blen_in: branch length info
    @param distn: state distribution at the root
    @param 
    @return: nx tree, map from edge to blen, map from edge to state
    """
    #XXX unfinished
    # populate the tree topology table
    edge_va_vb_list = (
            (0, 0, 5),
            (1, 1, 5),
            (2, 2, 6),
            (3, 3, 6),
            (4, 4, 7),
            (5, 6, 7),
            (6, 5, 7),
            )
    for edge_va_vb in edge_va_vb_list:
        cursor.execute('insert into topo values (?, ?, ?)', edge_va_vb)
    conn.commit()

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
    vertices = [t[0] for t in cursor]
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

    # build an undirected graph from the tree info
    G = nx.Graph()
    for edge, va, vb in edge_va_vb_list:
        G.add_edge(va, vb)

    # check that the graph is connected and has no cycles
    if not nx.is_connected(G):
        raise Exception('the tree is not connected')
    if nx.cycle_basis(G):
        raise Exception('the tree has a cycle')

    # build a directed breadth first tree starting at the distinguished vertex
    G_dag = nx.bfs_tree(G, root)

    # sample the unconditional histories
    conn = sqlite3.connect('histories.db')
    cursor = conn.cursor()
    cursor.execute(
            'create table histories ('
            'history integer, '
            'segment integer, '
            'va integer, '
            'vb integer, '
            'blen real, '
            'state integer, '
            'primary key (history, segment))')
    conn.commit()
    for i in range(nsamples):
        pass
        #XXX unfinished


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

