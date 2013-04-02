"""
This script is analogous to seq-gen.

But seq-gen is limited to nucleotide and amino acid alignments,
whereas we will allow more general alignments.
The continuous time Markov process is assumed to be time-reversible.
output schema
table alignment
offset integer, taxon integer, state integer,
primary key (offset, taxon)
"""

import argparse
import sqlite3

import numpy as np
import networkx as nx

import cmedbutil


#FIXME allow non-reversible rate matrices if the root of the tree is specified
#maybe use a different script for this.


#XXX copypasted
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

#XXX copypasted
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

#XXX copypasted and modified
def build_alignment_table(
        alignment_length, only_leaves,
        conn, root, G_dag, distn, states, rates, P):
    """
    Do unconditional forward sampling of an alignment at vertices of a tree.
    This is not conditional on the states at the leaves of the tree.
    It also does not use expm.
    This script throws out the transition histories along branches.
    @param alignment_length: sample this many columns of the alignment
    @param only_leaves: True if we only want to sample leaf states
    @param conn: database connection
    @param root: the arbitrary root of the directed acyclic graph
    @param G_dag: networkx directed acyclic graph
    @param distn: state distribution at the root
    @param states: ordered list of integer states
    @param rates: rates away from the states
    @param P: substitution distributions conditional on instantaneous change
    """

    # the table name is hardcoded
    table_name = 'alignment'

    # initialize the table
    cursor = conn.cursor()
    s = (
            'create table if not exists {table} ('
            'offset integer, '
            'taxon integer, '
            'state integer, '
            'primary key (offset, taxon))'
            ).format(table=table_name)
    cursor.execute(s)
    conn.commit()

    # get the list of taxa in the tree
    taxa = sorted(G_dag)

    # populate the table
    for offset in range(alignment_length):

        # sample the path history along the tree at an alignment column
        G_segmented = sample_history(root, G_dag, distn, rates, P)

        # insert the vertex states into the table
        for v in G_dag:

            # get edges adjacent to the vertex
            edges = [G_segmented[v][n] for n in G_segmented.neighbors(v)]

            # bail if we are at an internal vertex and we only want leaves
            if only_leaves and len(edges) > 1:
                continue

            # check that all of the edges have the same state
            uniqued_states = list(set(edge['state'] for edge in edges))
            if len(uniqued_states) != 1:
                raise Exception(
                    'expected all edges adjacent to a vertex '
                    'to share the same state')

            # extract the state shared by the edges adjacent to the vertex
            state = uniqued_states[0]
            
            # add to the table in the database
            s = 'insert into {table} values (?, ?, ?)'.format(table=table_name)
            t = (offset, v, state)
            cursor.execute(s, t)

        # commit the alignment column to the database
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

    #XXX backport this into the sample-unconditional-history script
    # read and validate the rate matrix info from the sqlite3 database file
    conn = sqlite3.connect(args.rates)
    cursor = conn.cursor()
    states, distn, Q = get_rate_matrix_info(cursor)
    conn.close()

    # Get a more convenient form of the rate matrix for forward simulation.
    rates, P = cmedbutil.decompose_rates(Q)

    #XXX backport this into the sample-unconditional-history script
    # extract the unrooted tree from the tree db file
    conn = sqlite3.connect(args.tree)
    cursor = conn.cursor()
    G = get_unrooted_tree(cursor)
    conn.close()

    # Pick the smallest vertex of G as an arbitrary root for sampling.
    # This choice can be made arbitrarily
    # because the rate matrix is currently defined to be time-reversible.
    root = min(G.nodes())

    # Build a directed breadth first tree starting at the distinguished vertex.
    # Note that the tree built by nx.bfs_tree and the edges yielded
    # by nx.bfs_edges do not retain the edge attributes.
    G_dag = nx.bfs_tree(G, root)
    for a, b in G_dag.edges():
        G_dag[a][b]['blen'] = G[a][b]['blen']

    # sample the unconditional columns of the alignment
    conn = sqlite3.connect(args.outfile)
    build_alignment_table(
            args.length, args.only_leaves,
            conn, root, G_dag, distn, states, rates, P)
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tree', default='brown.tree.db',
            help='unrooted tree as an sqlite3 database file')
    parser.add_argument('--rates', default='rate.matrix.db',
            help=('continuous-time Markov substitution process '
                'as an sqlite3 database file'))
    parser.add_argument('--length', type=cmedbutil.pos_int, default=10,
            help='sequence length')
    parser.add_argument('--only-leaves', action='store_true',
            help=('include only leaf states '
                'and do not include the states '
                'at the internal nodes of the tree'))
    parser.add_argument('-o', '--outfile', default='sampled.alignment.db',
            help='create this sqlite3 database file')
    main(parser.parse_args())

