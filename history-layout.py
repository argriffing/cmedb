"""
Associate 2d coordinates with each vertex of a sampled history.

This script should possibly be redesigned
for more suitable input and output.
"""

import argparse
import sqlite3
import math

import numpy as np
import networkx as nx


def main(args):

    # open the history db file
    # read the association of node index pairs to edge indices
    # read the association of edge indices to branch lengths
    conn = sqlite3.connect(args.infile)
    cursor = conn.cursor()
    cursor.execute('select va, vb, blen from history')
    va_vb_blen_list = list(cursor)
    conn.close()

    # I think that the equal arc layout does not depend on the
    # inital vertex chosen, and that this is some theorem.
    # Therefore we can start arbitrarily with the vertex with lowest index.
    va, vb, blen = zip(*va_vb_blen_list)
    vertices = sorted(set(va + vb))
    root = vertices[0]

    # build an undirected graph from the tree info
    G = nx.Graph()
    for va, vb, blen in va_vb_blen_list:
        G.add_edge(va, vb, blen=blen)

    # check that the graph is connected and has no cycles
    if not nx.is_connected(G):
        raise Exception('the tree is not connected')
    if nx.cycle_basis(G):
        raise Exception('the tree has a cycle')

    # Build a directed breadth first tree starting at the distinguished vertex.
    # Note that the tree built by nx.bfs_tree and the edges yielded
    # by nx.bfs_edges do not retain the edge attributes.
    G_dag = nx.bfs_tree(G, root)

    # For each vertex in the arbitrarily rooted tree,
    # compute the number of leaves in the subtree.
    vertex_to_subtree_nleaves = {}
    for node in reversed(nx.topological_sort(G_dag)):
        successors = G_dag.successors(node)
        if successors:
            nleaves = sum(vertex_to_subtree_nleaves[v] for v in successors)
        else:
            nleaves = 1
        vertex_to_subtree_nleaves[node] = nleaves

    # Compute the min and max angles allowed for the subtree at the root.
    # This allows roots of arbitrary degree, including 0 or 1.
    vertex_to_angles = {}
    root_successors = G_dag.successors(root)
    if len(root_successors) == 0:
        pass
    elif len(root_successors) == 1:
        nleaves = vertex_to_subtree_nleaves[root]
        vertex_to_angles[root] = np.array([1.0 / (nleaves + 1) , 2*math.pi])
    else:
        vertex_to_angles[root] = np.array([0, 2*math.pi])

    # For each non-root vertex in the arbitrarily rooted tree,
    # compute the min and max angles allowed for the subtree.
    for node in nx.topological_sort(G_dag):
        successors = G_dag.successors(node)
        if not successors:
            continue
        node_arc = vertex_to_angles[node]
        node_arcwidth = node_arc[1] - node_arc[0]
        node_nleaves = vertex_to_subtree_nleaves[node]
        low = node_arc[0]
        for v in successors:
            proportion = vertex_to_subtree_nleaves[v] / float(node_nleaves)
            high = low + proportion * node_arcwidth
            vertex_to_angles[v] = np.array([low, high])
            low = high

    # For each vertex in the arbitrarily rooted tree,
    # compute the x and y coordinates of the 2d layout.
    vertex_to_2d = {root : np.zeros(2)}
    for node in nx.topological_sort(G_dag):
        successors = G_dag.successors(node)
        for v in successors:
            theta = np.mean(vertex_to_angles[v])
            blen = G[node][v]['blen']
            d = np.array([np.cos(theta), np.sin(theta)])
            vertex_to_2d[v] = vertex_to_2d[node] + blen * d

    # open the database for writing
    # define the layout schema
    # populate table
    conn = sqlite3.connect(args.outfile)
    cursor = conn.cursor()
    cursor.execute(
            'create table layout ('
            'node integer, x real, y real, '
            'primary key (node))')
    conn.commit()
    for node in G:
        point = vertex_to_2d[node]
        x = float(point[0])
        y = float(point[1])
        triple = (node, x, y)
        cursor.execute('insert into layout values (?, ?, ?)', triple)
    conn.commit()
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', default='history.db',
            help='sampled state history as an sqlite3 database file')
    parser.add_argument('-o', '--outfile', default='layout.db',
            help='create this sqlite3 database file')
    main(parser.parse_args())

