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

import cmedbutil

def gen_branch_history_sample(state_in, blen_in, rates, P):
    """
    Path sampling along a branch with a known initial state.
    Yield (transition time, new state) pairs.
    The path sampling is conditional on the initial state
    but it is not conditional on the final state.
    So this is a 'forward' rather than a 'bridge' sampling.
    It does not require any matrix exponential computation.
    @param state_in: initial state
    @param blen_in: length of the branch
    @param rates: the rate away from each state
    @param P: transition matrix conditional on leaving a state
    """
    state = state_in
    nstates = len(rates)
    blen_accum = 0
    while True:
        rate = rates[state]
        scale = 1 / rate
        b = np.random.exponential(scale=scale)
        blen_accum += b
        if blen_accum >= blen_in:
            return
        distn = P[state]
        state = cmedbutil.random_category(distn)
        yield blen_accum, state


def sample_history(root, G_dag_in, distn, rates, P):
    """
    Sample state history on a rooted tree.
    The input tree is a networkx directed acyclic graph
    with 'blen' annotation on edges.
    @param root: root of the tree
    @param G_dag_in: tree as a networkx directed acyclic graph
    @param distn: state distribution at the root
    @param rates: the rate away from each state
    @param P: transition matrix conditional on leaving a state
    @return: networkx tree with edges annotated with 'blen' and 'state'
    """
    nstates = len(distn)
    if P.shape != (nstates, nstates):
        raise Exception('nstates mismatch')
    # Sample the initial state from the distribution.
    root_state = cmedbutil.random_category(distn)
    # Initialize the root state.
    vertex_to_state = {root : root_state}
    # Note that the subset of vertices that are shared with the
    # original unsegmented tree should have the same index
    # numbers in the new segmented tree.
    G = nx.Graph()
    vertices = list(G_dag_in)
    next_vertex = max(vertices) + 1
    for node in nx.topological_sort(G_dag_in):
        initial_state = vertex_to_state[node]
        for successor in G_dag_in.successors(node):
            blen = G_dag_in[node][successor]['blen']
            accum = 0
            segment_state = initial_state
            v = node
            for t, j in gen_branch_history_sample(root_state, blen, rates, P):
                G.add_edge(
                        v,
                        next_vertex,
                        blen=t-accum,
                        state=segment_state,
                        )
                segment_state = j
                v = next_vertex
                next_vertex += 1
                accum = t
            G.add_edge(
                    v,
                    successor,
                    blen=blen-accum,
                    state=segment_state,
                    )
            vertex_to_state[successor] = segment_state
    return G

def build_single_history_table(
        conn, table, root, G_dag, distn, states, rates, P):
    """
    @param conn: database connection
    @param table: validated alphanumeric table name
    @param root: root index of the tree
    @param G_dag: networkx directed acyclic graph
    @param distn: state distribution at the root
    @param states: ordered list of integer states
    @param rates: rates away from the states
    @param P: substitution distributions conditional on instantaneous change
    """
    # initialize the table
    cursor = conn.cursor()
    s = (
            'create table if not exists {table} ('
            'segment integer, '
            'va integer, '
            'vb integer, '
            'blen real, '
            'state integer, '
            'primary key (segment))'
            ).format(table=table)
    cursor.execute(s)
    conn.commit()
    # populate the table
    G_segmented = sample_history(root, G_dag, distn, rates, P)
    h_vertices = list(G_segmented)
    for segment_index, (va, vb) in enumerate(G_segmented.edges()):
        edge = G_segmented[va][vb]
        s = 'insert into %s values (?, ?, ?, ?, ?)' % table
        t = (
                segment_index,
                va,
                vb,
                edge['blen'],
                states[edge['state']],
                )
        cursor.execute(s, t)
    conn.commit()

def build_multiple_histories_table(
        nsamples, conn, table, root, G_dag, distn, states, rates, P):
    """
    @param nsamples: sample this many histories
    @param conn: database connection
    @param table: validated alphanumeric table name
    @param root: root index of the tree
    @param G_dag: networkx directed acyclic graph
    @param distn: state distribution at the root
    @param states: ordered list of integer states
    @param rates: rates away from the states
    @param P: substitution distributions conditional on instantaneous change
    """
    # initialize the table
    cursor = conn.cursor()
    s = (
            'create table if not exists {table} ('
            'history integer, '
            'segment integer, '
            'va integer, '
            'vb integer, '
            'blen real, '
            'state integer, '
            'primary key (history, segment))'
            ).format(table=table)
    cursor.execute(s)
    conn.commit()
    # populate the table
    for history_index in range(nsamples):
        G_segmented = sample_history(root, G_dag, distn, rates, P)
        h_vertices = list(G_segmented)
        for segment_index, (va, vb) in enumerate(G_segmented.edges()):
            edge = G_segmented[va][vb]
            s = 'insert into %s values (?, ?, ?, ?, ?, ?)' % table
            t = (
                    history_index,
                    segment_index,
                    va,
                    vb,
                    edge['blen'],
                    states[edge['state']],
                    )
            cursor.execute(s, t)
        conn.commit()

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



def main(args):

    # define the number of histories to sample
    nsamples = args.nsamples

    # validate table name
    if not args.table.isalnum():
        raise Exception('table name must be alphanumeric')

    # read and validate the rate matrix info from the sqlite3 database file
    conn = sqlite3.connect(args.rates)
    cursor = conn.cursor()
    states, distn, Q = get_rate_matrix_info(cursor)
    conn.close()

    # open the tree db file
    conn = sqlite3.connect(args.tree)
    cursor = conn.cursor()
    G = get_unrooted_tree(cursor)
    conn.close()

    # define the root node index for sampling;
    # the root is either user supplied,
    # or is taken to be the smallest node index otherwise
    vertices = sorted(G)
    if args.root is None:
        root = vertices[0]
    elif args.root in vertices:
        root = args.root
    else:
        raise Exception('the specified root should be a node index in the tree')

    # Build a directed breadth first tree starting at the distinguished vertex.
    # Note that the tree built by nx.bfs_tree and the edges yielded
    # by nx.bfs_edges do not retain the edge attributes.
    G_dag = nx.bfs_tree(G, root)
    for a, b in G_dag.edges():
        G_dag[a][b]['blen'] = G[a][b]['blen']

    # sample the unconditional history or histories
    rates, P = cmedbutil.decompose_rates(Q)
    conn = sqlite3.connect(args.outfile)
    if args.nsamples == 1:
        build_single_history_table(
                conn, args.table, root, G_dag, distn, states, rates, P)
    else:
        build_multiple_histories_table(
                nsamples,
                conn, args.table, root, G_dag, distn, states, rates, P)
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tree', default='brown.tree.db',
            help='unrooted tree as an sqlite3 database file')
    parser.add_argument('--rates', default='rate.matrix.db',
            help=('continuous-time Markov substitution process '
                'as an sqlite3 database file'))
    parser.add_argument('--root', type=int,
            help='root node index')
    parser.add_argument('--nsamples', type=cmedbutil.pos_int, default=5,
            help='sample this many histories')
    parser.add_argument('-o', '--outfile', default='sampled.histories.db',
            help='create this sqlite3 database file')
    parser.add_argument('--table', default='histories',
            help='name of table to create in the new database')
    main(parser.parse_args())

