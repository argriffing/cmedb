"""
Sample leaf-conditioned path histories from multiple alignment columns.

The output of this script is intended to be a huge file of sampled histories.
An arbitrary non-leaf root is picked for the tree.
For each aligned codon column in the alignment,
The subtree probability is computed for each node using dynamic programming,
from the leaves of the tree up to the root.
Then state distributions are computed for each node,
for each possible parent node state.
These conditional distributions are used to jointly sample
states at the internal nodes of the tree,
starting from the arbitrary root node.
Then endpoint-conditioned path histories are sampled along each branch.
.
Because this script needs expm to sample states at internal vertices,
we need to map between nominal state indices and matrix state indices.
.
For now, do not treat nhistories=1 as a special case.
"""


import argparse
import sqlite3
from collections import defaultdict

import numpy as np
import scipy.linalg

import cmedbutil


#XXX copypasted
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
    cmedbutil.assert_stochastic_vector(distn)

    # assert that the rate matrix is actually a rate matrix
    cmedbutil.assert_rate_matrix(Q)

    # assert that the distribution is at equilibrium w.r.t. the rate matrix
    cmedbutil.assert_equilibrium(Q, distn)

    # assert that the detailed balance equations are met
    cmedbutil.assert_detailed_balance(Q, distn)

    # return the validated inputs describing the stochastic process
    return states, distn, Q


#XXX copypasted
def get_unrooted_tree(cursor):
    """
    Extract an unrooted tree from an sqlite3 database connection.
    @return: undirected networkx graph with blen edge attributes
    """

    # read the association of node index pairs to edge indices
    # read the association of edge indices to branch lengths
    cursor.execute('select edge, va, vb from topo')
    edge_va_vb_list = sorted(cursor)
    cursor.execute('select edge, blen from blen')
    edge_to_blen = dict(cursor)

    # build an undirected graph from the tree info
    G = nx.Graph()
    for edge, va, vb in edge_va_vb_list:
        G.add_edge(va, vb, blen=edge_to_blen[edge])

    # check that the graph is connected and has no cycles
    cmedbutil.assert_connected_acyclic_graph(G)

    # return the graph
    return G


def sample_site_history(foo, bar):
    """
    Endpoint conditioned path history sampling at a site.
    All inputs and outputs of this function should
    represent genetic states using row indices into a matrix,
    as opposed to the ascii name of the state
    or the numerical identifier used in the sqlite
    rate matrix and alignment files.
    """
    # Use a Felsenstein-like algorithm
    # to compute subtree probs given a subtree root state.
    # This depends on the rooted tree structure, the edge_to_P map,
    # and the states of the alignment column at the leaves.
    v_to_subtree_probs = {}


def main(args):

    # define the number of histories to sample
    nsamples = args.nsamples

    # read and validate the rate matrix info from the sqlite3 database file
    conn = sqlite3.connect(args.rates)
    cursor = conn.cursor()
    states, distn, Q = get_rate_matrix_info(cursor)
    conn.close()

    # read the unrooted tree from the database
    conn = sqlite3.connect(args.tree)
    cursor = conn.cursor()
    G = get_unrooted_tree(cursor)
    conn.close()

    # read the alignment from the database
    conn = sqlite3.connect(args.alignment)
    cursor = conn.cursor()
    cursor.execute(
            'select offset, taxon, state from alignment '
            'order by offset, taxon')
    offset_taxon_state = list(cursor)
    cursor.execute('select offset from alignment')
    offsets = set(t[0] for t in cursor)
    conn.close()

    # define the map from state to rate matrix index
    s_to_i = dict((s, i) for i, s in enumerate(states))

    # for each offset get a map from leaf taxon to genetic state
    offset_to_leaf_map = defaultdict(dict)
    for offset, taxon, state in offset_taxon_state:
        offset_to_leaf_map[offset][taxon] = state

    # Define the leaf and non-leaf vertices.
    # Pick an arbitrary non-leaf vertex to act as the root.
    leaves = sorted(v for v in G if G.degree(v) == 1)
    non_leaves = sorted(set(G) - set(leaves))
    root = non_leaves[0]

    # Build a directed breadth first tree starting at the distinguished vertex.
    # Note that the tree built by nx.bfs_tree and the edges yielded
    # by nx.bfs_edges do not retain the edge attributes.
    G_dag = nx.bfs_tree(G, root)
    for a, b in G_dag.edges():
        G_dag[a][b]['blen'] = G[a][b]['blen']

    # Use expm to associate a conditional transition matrix
    # with each directed branch in the tree.
    # This needs to be done only once for the whole simulation,
    # regardless of the number of requested histories
    # and regardless of the number of columns in the alignment.
    edge_to_P = {}
    for a, b in G_dag.edges():
        edge_to_P[a, b] = scipy.linalg.expm(Q * G_dag[a][b]['blen'])

    # Get a decomposition of the rate matrix for path sampling on branches.
    rates, P = cmedbutil.decompose_rates(Q)

    # open a database file for writing
    conn = sqlite3.connect(args.outfile)
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
            'primary key (history, offset, segment))'
            ).format(table=table)
    cursor.execute(s)

    # populate the database
    for history_index in range(args.nhistories):
        for offset in offsets:
            for seg, va, vb, blen, state in sample_site_history(
                    foo, bar):
                s = 'insert into histories values (?, ?, ?, ?, ?, ?, ?)'
                t = (history_index, offset, seg, va, vb, blen, state)
                cursor.execute(s, t)
                conn.commit()



    # close the output database
    conn.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tree', default='brown.tree.db',
            help='unrooted tree as an sqlite3 database file')
    parser.add_argument('--rates', default='rate.matrix.db',
            help=('reversible continuous-time Markov substitution process '
                'as an sqlite3 database file'))
    parser.add_argument('--alignment', default='alignment.db',
            help='input i.i.d. aligned leaf states as an sqlite3 database file')
    parser.add_argument('--nhistories', type=cmedbutil.pos_int, default=5,
            help='sample this many histories')
    parser.add_argument('--outfile', default='sampled.histories.db',
            help='create this sqlite3 database file')
    main(parser.parse_args())

