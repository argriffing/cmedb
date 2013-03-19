"""
Identify the isostate subtrees within a sampled substitution history.

The input is an unrooted tree for which each edge has been associated
with a state.
The output associated two more labels with each edge.
The first label is an isostate label.
Two edges in the tree have the same isostate label
if and only if there exists a path between the edges
such that the state does not change along the path.
The second label is 0 or 1.
It is actually associated with the isostate label.
schema
.
table isostate
va integer, vb integer, isostate integer
.
table bipartite
isostate integer, parity integer
"""

import argparse
import sqlite3

import networkx as nx


def get_state_segmentation(G_in):
    """
    Segment the tree according to state.
    This does not use the branch lengths or the layout.
    @param G_in: undirected graph with state annotation on edges
    @return: va_vb_isostate_list, isostate_to_parity
    """

    # get leaf vertices
    # pick an arbitrary leaf as a distinguished (root) vertex
    vertices = list(G_in)
    leaves = sorted(v for v in vertices if len(G_in.neighbors(v)) == 1)
    root = leaves[0]

    # Build a directed breadth first tree starting at the distinguished vertex.
    # Note that the tree built by nx.bfs_tree and the edges yielded
    # by nx.bfs_edges do not retain the edge attributes.
    G_dag = nx.bfs_tree(G_in, root)

    # initialize the tree of isostate adjacencies
    G_isostate = nx.Graph()

    # Each contig is defined by a set of edges.
    root_state = G_in[root][G_dag.successors(root)[0]]['state']
    root_edge_list = []
    root_contig_index = 0
    contig_states = [root_state]
    contig_edge_lists = [root_edge_list]
    vertex_to_contig_index = {root : root_contig_index}
    for node in nx.topological_sort(G_dag):
        ci = vertex_to_contig_index[node]
        ci_state = contig_states[ci]
        successors = G_dag.successors(node)
        for v in successors:
            state = G_in[node][v]['state']
            if state == ci_state:
                contig_edge_lists[ci].append((node, v))
                vertex_to_contig_index[v] = ci
            else:
                ci_next = len(contig_states)
                G_isostate.add_edge(ci, ci_next)
                vertex_to_contig_index[v] = ci_next
                contig_states.append(state)
                contig_edge_lists.append([(node, v)])

    # Convert the G_isostate graph into a map from
    # isostate labels to parities.
    isostate_to_parity = {0 : 0}
    for va, vb in nx.bfs_edges(G_isostate, 0):
        isostate_to_parity[vb] = 1 - isostate_to_parity[va]
    
    # Get the isostate label associated with each edge.
    va_vb_isostate_list = []
    for isostate_label, edge_list in enumerate(contig_edge_lists):
        for va, vb in edge_list:
            va_vb_isostate_list.append((va, vb, isostate_label))
    return va_vb_isostate_list, isostate_to_parity


def main(args):

    # open the history db file
    # read the segmented tree graph
    conn = sqlite3.connect(args.history)
    cursor = conn.cursor()
    cursor.execute('select va, vb, state from history')
    va_vb_state_list = list(cursor)
    conn.close()

    # construct a list of vertices
    va_list, vb_list, state_list = zip(*va_vb_state_list)
    vertices = sorted(set(va_list + vb_list))

    # some input validation
    if len(vertices) < 2:
        raise Exception('expected at least one edge')

    # build an undirected graph from the tree info
    G = nx.Graph()
    for va, vb, state in va_vb_state_list:
        G.add_edge(va, vb, state=state)

    # check that the graph is connected and has no cycles
    if not nx.is_connected(G):
        raise Exception('the tree is not connected')
    if nx.cycle_basis(G):
        raise Exception('the tree has a cycle')

    # Check that all edges sharing a high degree vertex
    # have the same state.
    for node in vertices:
        neighbors = G.neighbors(node)
        if len(neighbors) >= 3:
            state_set = set(G[node][v]['state'] for v in neighbors)
            if len(state_set) > 1:
                raise Exception(
                        'expected all edges adjacent to a high degree node '
                        'to share the same state')

    # get the segmentation
    va_vb_isostate_list, isostate_to_parity = get_state_segmentation(G)

    # open the database for writing
    # define the layout schema
    # populate table
    conn = sqlite3.connect(args.segmentation)
    cursor = conn.cursor()
    cursor.execute(
            'create table isostate ('
            'va integer, vb integer, isostate integer, '
            'primary key (va, vb))')
    cursor.execute(
            'create table bipartite ('
            'isostate integer, parity integer, '
            'primary key (isostate))')
    conn.commit()
    for triple in va_vb_isostate_list:
        cursor.execute('insert into isostate values (?, ?, ?)', triple)
    conn.commit()
    for pair in isostate_to_parity.items():
        cursor.execute('insert into bipartite values (?, ?)', pair)
    conn.commit()
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--history', default='history.db',
            help='input sampled state history as an sqlite3 database file')
    parser.add_argument('--segmentation', default='segmentation.db',
            help='output segmentation as an sqlite3 database file')
    main(parser.parse_args())

