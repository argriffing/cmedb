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
"""

#FIXME redo this to not use numpy,
# and to more flexibly support sparse rate matrices
# without worrying about extra conversions between states
# and indices into ndarrays.
# This will mean not using the usual rate matrix reader.

import argparse
import sqlite3

import numpy as np

import cmedbutil


# XXX under construction
def gen_branch_history_sample(
        primary_state_in, blink_states_in, blen_in,
        primary_Q,
        partition, on_rate, off_rate,
        ):
    """
    This function is defined in analogy to gen_branch_history_sample.
    Path sampling along a branch with a known initial state.
    Yield (transition time, new state) pairs.
    This function takes arguments divided into three groups.
    The first group defines the initial state and the path duration.
    The second group defines the primary process
    without regard to blinking.
    The third group defines the partition of the primary state space
    and the toggling rates of the modulating blinking processes.
    @param primary_state_in: initial integer state of the observable process
    @param blink_state_in: initial binary tuple state of the blinking process
    @param blen_in: length of the branch
    @param primary_Q: primary process rate matrix without regard to blinking
    @param partition: maps a primary state to a partition index for blinking
    @param on_rate: instantaneous off-to-on blinking rate
    @param off_rate: instantaneous on-to-off blinking rate
    """
    #XXX the docstring is up to date but the implmenetation is not
    primary_state = primary_state_in
    blink_state = tuple(blink_state_in)
    blen_accum = 0
    while True:

        # Compute the total rate out of the compound state.
        # This is the sum of allowed primary transition rates
        # and allowed blink toggle rates.
        total_primary_rate = 0.0
        for j in primary_Q:
            pass

        rate = rates[state]
        scale = 1 / rate
        b = np.random.exponential(scale=scale)
        blen_accum += b
        if blen_accum >= blen_in:
            return
        distn = P[state]
        state = cmedbutil.random_category(distn)
        yield blen_accum, state


# XXX this needs to be changed or deleted
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

# XXX this might need to be changed
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


# XXX this has been copypasted and needs to be changed
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


#XXX the command line args should be correct now
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

