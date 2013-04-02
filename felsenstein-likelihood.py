"""
Compute a likelihood using the Felsenstein pruning algorithm.

This likelihood is computed given an unrooted phylogenetic tree
with branch lengths,
aligned columns of states at the leaves of the tree,
and a rate matrix defining a continuous time-reversible Markov process.
output is just a likelhood written to the screen
or a log likelihood or something.
maybe an average log likelhood.
"""

import argparse
import sqlite3
from collections import defaultdict

import numpy as np
import networkx as nx

import cmedbutil

#XXX XXX this script will depend on an alignment -> patterns preprocessing step


#TODO allow a root and a non-reversible rate matrix

#TODO use fancy things like the BEAGLE library or SSE or GPU



#XXX copypasted from slowedml/alignll
def _ll_helper(
        ov, v_to_children, de_to_P, root_prior,
        patterns, pat_mults,
        fn,
        ):
    """
    This is purely a helper function.
    It should not be called outside of its own module.
    The P matrices and the root prior may be algopy objects.
    @param ov: ordered vertices with child vertices before parent vertices
    @param v_to_children: map from a vertex to a sequence of child vertices
    @param de_to_P: map from a directed edge to a transition matrix
    @param root_prior: equilibrium distribution at the root
    @param patterns: each pattern assigns a state to each leaf
    @param pat_mults: a multiplicity for each pattern
    @param fn: a per-pattern log likelihood evaluation function
    @return: log likelihood
    """
    npatterns = patterns.shape[0]
    lls = algopy.zeros(
            npatterns,
            dtype=de_to_P.values()[0],
            )
    for i in range(npatterns):
        lls[i] = fn(ov, v_to_children, patterns[i], de_to_P, root_prior)
    return algopy.dot(lls, pat_mults)


#XXX copypasted from slowedml/sitell
def brute_log_likelihood(ov, v_to_children, pattern, de_to_P, root_prior):
    """
    Brute force likelihood calculation.
    @param ov: ordered vertices with child vertices before parent vertices
    @param v_to_children: map from a vertex to a sequence of child vertices
    @param pattern: an array that maps vertex to state, or to -1 if internal
    @param de_to_P: map from a directed edge to a transition matrix
    @param root_prior: equilibrium distribution at the root
    @return: log likelihood
    """
    nvertices = len(pattern)
    nstates = len(root_prior)
    root = ov[-1]
    v_unknowns = [v for v, state in enumerate(pattern) if state == -1]
    n_unknowns = len(v_unknowns)

    # Construct the set of directed edges on the tree.
    des = set((p, c) for p, cs in v_to_children.items() for c in cs)

    # Compute the likelihood by directly summing over all possibilities.
    likelihood = 0
    for assignment in itertools.product(range(nstates), repeat=n_unknowns):

        # Fill in the state assignments for all vertices.
        augmented_pattern = np.array(pattern)
        for v, state in zip(v_unknowns, assignment):
            augmented_pattern[v] = state

        # Add to the log likelihood.
        edge_prob = 1.0
        for p, c in des:
            p_state = augmented_pattern[p]
            c_state = augmented_pattern[c]
            edge_prob *= de_to_P[p, c][p_state, c_state]
        likelihood += root_prior[augmented_pattern[root]] * edge_prob

    # Return the log likelihood.
    return algopy.log(likelihood)


#XXX copypasted from slowedml/sitell
def felsenstein_log_likelihood(ov, v_to_children, pattern, de_to_P, root_prior):
    """
    Felsenstein pruning likelihood calculation using dynamic programming.
    @param ov: ordered vertices with child vertices before parent vertices
    @param v_to_children: map from a vertex to a sequence of child vertices
    @param pattern: an array that maps vertex to state, or to -1 if internal
    @param de_to_P: map from a directed edge to a transition matrix
    @param root_prior: equilibrium distribution at the root
    @return: log likelihood
    """
    nvertices = len(ov)
    nstates = len(root_prior)
    states = range(nstates)
    root = ov[-1]

    # Initialize the map from vertices to subtree likelihoods.
    likelihoods = np.ones(
            (nvertices, nstates),
            dtype=de_to_P.values()[0],
            )

    # Compute the subtree likelihoods using dynamic programming.
    for v in ov:
        for pstate in range(nstates):
            for c in v_to_children.get(v, []):
                P = de_to_P[v, c]
                likelihoods[v, pstate] *= np.dot(P[pstate], likelihoods[c])
        state = pattern[v]
        if state >= 0:
            for s in range(nstates):
                if s != state:
                    likelihoods[v, s] = 0

    # Get the log likelihood by summing over equilibrium states at the root.
    return np.log(np.dot(root_prior, likelihoods[root]))



def get_leaf_alignment(cursor):
    """
    Extract a leaf state alignment from an sqlite3 database cursor.
    @param cursor: database cursor
    @return: sequence of maps from taxon to state
    """

    # read the data from the database cursor
    cursor.execute('select offset, taxon, state from alignment')
    offset_to_pairs = defaultdict(list)
    for offset, taxon, state in cursor:
        offset_to_pairs[offset].append((taxon, state))
    assignments = [dict(pairs) for pairs in offset_to_pairs.values()]
    return assignments


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


def main(args):

    #XXX backport this into the sample-unconditional-history script
    # read and validate the rate matrix info from the sqlite3 database file
    conn = sqlite3.connect(args.rates)
    cursor = conn.cursor()
    states, distn, Q = get_rate_matrix_info(cursor)
    conn.close()

    #XXX backport this into the sample-unconditional-history script
    # extract the unrooted tree from the tree db file
    conn = sqlite3.connect(args.tree)
    cursor = conn.cursor()
    G = get_unrooted_tree(cursor)
    conn.close()

    # Extract the leaf nodes from the undirected tree.
    leaves = sorted(v for v in G if G.degree(v) < 2)
    if not leaves:
        raise Exception('could not find any leaves')

    # To be unnecessarily restrictive,
    # we will require that the root is an internal vertex in the tree.
    # The choice can be made arbitrarily
    # because the rate matrix is currently defined to be time-reversible.
    best_degree, best_vertex = max((G.degree(v), v) for v in G)
    root = best_vertex

    # Build a directed breadth first tree starting at the distinguished vertex.
    # Note that the tree built by nx.bfs_tree and the edges yielded
    # by nx.bfs_edges do not retain the edge attributes.
    G_dag = nx.bfs_tree(G, root)
    for a, b in G_dag.edges():
        G_dag[a][b]['blen'] = G[a][b]['blen']

    # Read the leaf alignment from the database.
    conn = sqlite3.connect(args.leaf_alignment)
    cursor = conn.cursor()
    assignments = get_leaf_alignment(cursor)
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tree', default='brown.tree.db',
            help='unrooted tree as an sqlite3 database file')
    parser.add_argument('--rates', default='rate.matrix.db',
            help=('continuous-time Markov substitution process '
                'as an sqlite3 database file'))
    parser.add_argument('--leaf-alignment', default='leaf.alignment.db',
            help='aligned states at the leaves of the tree')
    main(parser.parse_args())

