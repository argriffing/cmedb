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

def decompose_rates(Q):
    """
    Break a rate matrix into two parts.
    The first part consists of the rates away from each state;
    this information is contained in the diagonal of the rate matrix.
    The second part consists of a transition matrix
    that defines the distribution over sink states conditional
    on an instantaneous change away from a given source state.
    Note that this function never requires expm of Q.
    Also, P preserves the sparsity pattern of Q.
    @param Q: rate matrix
    @return: rates, P
    """
    nstates = len(Q)
    rates = -np.diag(Q)
    P = np.array(Q)
    for i, rate in enumerate(rates):
        if rate:
            P[i, i] = 0
            P[i] /= rate
    return rates, P

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
        new_state = np.dot(np.arange(nstates), np.random.multinomial(1, distn))
        yield blen_accum, new_state


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
    # Python and numpy and scipy seem to not have categorial distribution
    # (weighted choice) implementations, so I will use the multinomial
    # sampling function instead.
    v = np.random.multinomial(1, distn)
    root_state = np.dot(np.arange(nstates), v)
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
    # construct the distribution of states at the root
    cursor.execute('select state, prob from distn')
    state_prob_pairs = list(cursor)
    distn = np.zeros(nstates)
    for state, prob in state_prob_pairs:
        distn[s_to_i[state]] = prob
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
        G.add_edge(va, vb, blen=edge_to_blen[edge])

    # check that the graph is connected and has no cycles
    if not nx.is_connected(G):
        raise Exception('the tree is not connected')
    if nx.cycle_basis(G):
        raise Exception('the tree has a cycle')

    # Build a directed breadth first tree starting at the distinguished vertex.
    # Note that the tree built by nx.bfs_tree and the edges yielded
    # by nx.bfs_edges do not retain the edge attributes.
    G_dag = nx.bfs_tree(G, root)
    for a, b in G_dag.edges():
        G_dag[a][b]['blen'] = G[a][b]['blen']

    # sample the unconditional histories
    conn = sqlite3.connect('histories.db')
    cursor = conn.cursor()
    cursor.execute(
            'create table if not exists histories ('
            'history integer, '
            'segment integer, '
            'va integer, '
            'vb integer, '
            'blen real, '
            'state integer, '
            'primary key (history, segment))')
    conn.commit()

    # get the samples
    print 'histories:'
    rates, P = decompose_rates(Q)
    for i in range(nsamples):
        G_segmented = sample_history(root, G_dag, distn, rates, P)
        print nx.to_edgelist(G_segmented)
        print


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tree', required=True,
            help='unrooted tree as an sqlite3 database file')
    parser.add_argument('--rates', required=True,
            help=('continuous-time Markov substitution process '
                'as an sqlite3 database file'))
    parser.add_argument('--root', type=int,
            help='root node index')
    parser.add_argument('--nsamples', type=int, default=5,
            help='sample this many histories')
    main(parser.parse_args())

