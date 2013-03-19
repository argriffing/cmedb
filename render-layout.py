"""
Create an image file visualization of a state history on a tree using pycairo.

The substitution history and the layout should be provided
as sqlite3 database files.
"""

import argparse
import sqlite3
import cairo

import numpy as np
import networkx as nx

def get_state_segmentation(va_vb_state_list):
    """
    Segment the tree according to state.
    This does not use the branch lengths or the layout.
    """

    # define the vertices and an arbitrary root
    va_list, vb_list, state_list = zip(*va_vb_state_list)
    vertices = sorted(set(va_list + vb_list))
    root = vertices[0]
    
    # construct the unrooted networkx tree
    G = nx.Graph()
    for va, vb, state in va_vb_state_list:
        G.add_edge(va, vb, state=state)

    # check that the graph is connected and has no cycles
    if not nx.is_connected(G):
        raise Exception('the tree is not connected')
    if nx.cycle_basis(G):
        raise Exception('the tree has a cycle')

    # Build a directed breadth first tree starting at the distinguished vertex.
    # Note that the tree built by nx.bfs_tree and the edges yielded
    # by nx.bfs_edges do not retain the edge attributes.
    G_dag = nx.bfs_tree(G, root)

    # Each contig is defined by a set of edges and a parity.
    root_state = G[root][G_dag.successors(root)[0]]['state']
    root_parity = 0
    root_edge_list = []
    root_contig_index = 0
    contig_states = [root_state]
    contig_parities = [root_parity]
    contig_edge_lists = [root_edge_list]
    vertex_to_contig_index = {root : root_contig_index}
    for node in nx.topological_sort(G_dag):
        ci = vertex_to_contig_index[node]
        ci_state = contig_states[ci]
        successors = G_dag.successors(node)
        for v in successors:
            state = G[node][v]['state']
            if state == ci_state:
                contig_edge_lists[ci].append((node, v))
                vertex_to_contig_index[v] = ci
            else:
                vertex_to_contig_index[v] = len(contig_states)
                contig_states.append(state)
                contig_parities.append(1 - contig_parities[ci])
                contig_edge_lists.append([(node, v)])
    
    # For the purposes of rendering,
    # we only need to annotate edges with the parity.
    va_vb_parity_list = []
    for parity, edge_list in zip(contig_parities, contig_edge_lists):
        for va, vb in edge_list:
            va_vb_parity_list.append((va, vb, parity))
    return va_vb_parity_list


def main(args):

    # open the history db file
    # read the segmented tree graph
    conn = sqlite3.connect(args.history)
    cursor = conn.cursor()
    cursor.execute('select va, vb, blen, state from history')
    va_vb_blen_state_list = list(cursor)
    conn.close()

    # open the layout db file
    # read the 2d layout of the vertices
    conn = sqlite3.connect(args.layout)
    cursor = conn.cursor()
    cursor.execute('select node, x, y from layout')
    node_x_y_list = list(cursor)
    conn.close()

    # get a parity associated with each contiguous state
    va_vb_state_list = [
            (va, vb, state) for va, vb, blen, state in va_vb_blen_state_list]
    va_vb_parity_list = get_state_segmentation(va_vb_state_list)

    # define the 2d location of each vertex
    vertex_to_2d = dict((v, np.array([x, y])) for v, x, y in node_x_y_list)

    # define width and height of screen
    screenwidth = 400.0
    screenheight = 200.0

    # get bounding box
    xmin = min(p[0] for p in vertex_to_2d.values())
    ymin = min(p[1] for p in vertex_to_2d.values())
    xmax = max(p[0] for p in vertex_to_2d.values())
    ymax = max(p[1] for p in vertex_to_2d.values())

    # get center of bounding box
    xmid = (xmin + xmax) / 2.0
    ymid = (ymin + ymax) / 2.0
    mid = np.array([xmid, ymid])

    # rescale by this much
    xscale = screenwidth / (xmax - xmin)
    yscale = screenheight / (ymax - ymin)
    scale = min(xscale, yscale)

    # Construct an undirected tree graph
    # and decorate its edges with attributes.
    G = nx.Graph()
    for va, vb, blen, state in va_vb_blen_state_list:
        G.add_edge(va, vb)
        G[va][vb]['blen'] = blen
        G[va][vb]['state'] = state
    for va, vb, parity in va_vb_parity_list:
        print va, vb, parity
        G[va][vb]['parity'] = parity

    # draw the image
    surface = cairo.SVGSurface('out.svg', screenwidth, screenheight)
    ctx = cairo.Context(surface)
    ctx.translate(screenwidth/2, screenheight/2)
    #for va, vb, blen, state in va_vb_blen_state_list:
    for va, vb in G.edges():
        print va, vb
        parity = G[va][vb]['parity']
        pa = (vertex_to_2d[va] - mid) * scale
        pb = (vertex_to_2d[vb] - mid) * scale
        if parity:
            ctx.set_source_rgb(1.0, 0.6, 0.6)
        else:
            ctx.set_source_rgb(0.6, 0.6, 1.0)
        ctx.move_to(pa[0], pa[1])
        ctx.line_to(pb[0], pb[1])
        ctx.stroke()
    ctx = None
    surface.finish()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--history', default='history.db',
            help='sampled state history as an sqlite3 database file')
    parser.add_argument('--layout', default='layout.db',
            help='2d node layout as an sqlite3 database file')
    #parser.add_argument('--imageformat', choices=('png', 'svg', 'pdf'),
            #help='choose among image formats')
    parser.add_argument('-o', '--outfile', default='out.png',
            help='create this image file')
    main(parser.parse_args())

