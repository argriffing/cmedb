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
import random
import math
from collections import defaultdict

import numpy as np
import networkx as nx
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


# XXX similar to code in felsenstein-likelihood.py
def get_v_to_subtree_probs(nstates, G_dag, edge_to_P, leaf_to_state_index):
    """
    Get partial likelihoods on an arbitrarily rooted tree.
    Use a Felsenstein-like algorithm
    to compute subtree probs given a subtree root state.
    This depends on the rooted tree structure, the edge_to_P map,
    and the states of the alignment column at the leaves.
    @param nstates: number of states
    @param G_dag: arbitrarily rooted genetic history with branch lengths
    @param edge_to_P: (a, b) to transition matrix
    @param leaf_to_state_index: alignment column information
    @return: map from a vertex to a vector of state-conditioned subtree probs
    """
    v_to_subtree_probs = {}
    for a in reversed(nx.topological_sort(G_dag)):
        subtree_probs = np.ones(nstates, dtype=float)
        successors = G_dag.successors(a)
        if successors:
            for state_index in range(nstates):
                for b in successors:
                    P = edge_to_P[a, b]
                    p = np.dot(P[state_index], v_to_subtree_probs[b])
                    subtree_probs[state_index] *= p
        else:
            state_index = leaf_to_state_index[a]
            subtree_probs[state_index] = 1
        v_to_subtree_probs[a] = subtree_probs
    return v_to_subtree_probs


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


#XXX copypasted
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

#XXX copypasted
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


def sample_site_history(Q, distn, G_dag, edge_to_P, leaf_to_state_index):
    """
    Endpoint conditioned path history sampling at a site.
    All inputs and outputs of this function should
    represent genetic states using row indices into a matrix,
    as opposed to the ascii name of the state
    or the numerical identifier used in the sqlite
    rate matrix and alignment files.
    @param Q: rate matrix
    @param distn: equilibrium distribution of the reversible process
    @param G_dag: arbitrarily rooted genetic history with branch lengths
    @param edge_to_P: (a, b) to transition matrix
    @param leaf_to_state_index: alignment column information
    @return: path history as a list of tuples
    """

    # Define the leaf and non-leaf vertices.
    # Pick an arbitrary non-leaf vertex to act as the root.
    leaves = sorted(v for v in G_dag if G_dag.degree(v) == 1)
    non_leaves = sorted(set(G_dag) - set(leaves))
    root = non_leaves[0]
    
    # Define more temporary info.
    leaf_set = set(leaves)
    nleaves = len(leaf_set)
    nstates = Q.shape[0]

    # Use some dynamic programming to get the
    # the subtree probability for each internal node
    # for each state index at that internal node.
    v_to_subtree_probs = get_v_to_subtree_probs(
            nstates, G_dag, edge_to_P, leaf_to_state_index)

    # Sample the states at internal vertices,
    # conditional on the values at the leaves
    # and on the parent state.
    v_to_state_index = {}
    for v in nx.topological_sort(G_dag):
        predecessors = G_dag.predecessors(v)
        successors = G_dag.successors(v)
        if successors:
            if not predecessors:
                prior_distn = distn
            elif len(predecessors) == 1:
                parent = predecessors[0]
                parent_state_index = v_to_state_index[parent]
                P = edge_to_P[parent, v]
                prior_distn = P[parent_state_index]
            else:
                raise Exception('found a node with multiple predecessors')
            weights = prior_distn * v_to_subtree_probs[v]
            category_distn = weights / np.sum(weights)
            v_to_state_index[v] = cmedbutil.random_category(category_distn)
        else:
            v_to_state_index[v] = leaf_to_state_index[v]

    # Get a decomposition of the rate matrix for path sampling on branches.
    rates, P = cmedbutil.decompose_rates(Q)

    # Endpoint conditioned path sampling along each branch.
    # Create new vertices at substitution points along the paths.
    # Each segment gets its own index.
    site_history = []
    next_vertex = max(G_dag) + 1
    next_seg = 0
    for a, b in G_dag.edges():
        initial_si = v_to_state_index[a]
        final_si = v_to_state_index[b]
        total_length = G_dag[a][b]['blen']
        si_blen_pairs = get_modified_rejection_sample(
                total_length, initial_si, final_si, rates, P)
        nsegs = len(state_blen_pairs)
        for i in range(nsegs-1):
            path_vertices.append(next_vertex)
            next_vertex += 1
        path_vertices == [a] + path_vertices + [b]
        path_edges = list(cmedbutil.pairwise(path_vertices))
        for (va, vb), (si, blen) in zip(path_edges, si_blen_pairs):
            t = (next_seg, va, vb, blen, si)
            site_history.append(t)
            next_seg += 1

    # return the sampled history
    return site_history


def main(args):

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
            'primary key (history, offset, segment))')
    cursor.execute(s)

    # populate the database
    for history_index in range(args.nhistories):
        for offset in offsets:
            if args.verbose:
                print 'sampling a new (history, offset)'
                print 'history index:', history_index
                print 'offset:', offset
                print
            leaf_to_state_index = {}
            for leaf in leaves:
                state = offset_to_leaf_map[offset][leaf]
                leaf_to_state_index[leaf] = s_to_i[state]
            for seg, va, vb, blen, state_index in sample_site_history(
                    Q, distn, G_dag, edge_to_P, leaf_to_state_index):
                state = states[state_index]
                s = 'insert into histories values (?, ?, ?, ?, ?, ?, ?)'
                t = (history_index, offset, seg, va, vb, blen, state)
                cursor.execute(s, t)
                conn.commit()

    # close the output database
    conn.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--verbose', action='store_true',
            help='spam more text')
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

