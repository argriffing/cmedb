"""
Create an image file visualization of a state history on a tree using pycairo.

The substitution history and the layout should be provided
as sqlite3 database files.
"""

import argparse
import sqlite3
from collections import defaultdict

import cairo
import numpy as np
import networkx as nx


def main(args):

    # open the history db file
    # read the segmented tree graph
    conn = sqlite3.connect(args.history)
    cursor = conn.cursor()
    seg_va_vb = list(cursor.execute('select segment, va, vb from history'))
    seg_to_blen = dict(cursor.execute('select segment, blen from history'))
    seg_to_state = dict(cursor.execute('select segment, state from history'))
    conn.close()

    # open the segmentation db file
    # read some more annotations about the segmentation
    conn = sqlite3.connect(args.segmentation)
    cursor = conn.cursor()
    cursor.execute('select segment, isostate from isostate')
    seg_to_isostate = dict(cursor)
    cursor.execute('select isostate, parity from bipartite')
    isostate_to_parity = dict(cursor)
    conn.close()

    # open the layout db file
    # read the 2d layout of the vertices
    conn = sqlite3.connect(args.layout)
    cursor = conn.cursor()
    cursor.execute('select node, x, y from layout')
    node_x_y_list = list(cursor)
    conn.close()

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
    for seg, va, vb in seg_va_vb:
        blen = seg_to_blen[seg]
        state = seg_to_state[seg]
        parity = isostate_to_parity[seg_to_isostate[seg]]
        G.add_edge(va, vb, segment=seg, blen=blen, state=state, parity=parity)

    # define colors associated with edges
    black_rgb = (0.0, 0.0, 0.0)
    light_red_rgb = (1.0, 0.6, 0.6)
    light_blue_rgb = (0.6, 0.6, 1.0)
    parity_to_rgb = {0 : light_red_rgb, 1 : light_blue_rgb}

    # Get the text labels associated with the edges.
    # Most edges will have no associated text label.
    # The longest segment of each isostate subtree
    # will have a text label that corresponds to the state.
    isostate_to_blen_seg_pairs = defaultdict(list)
    for seg, blen in seg_to_blen.items():
        isostate = seg_to_isostate[seg]
        blen_seg_pair = (blen, seg)
        isostate_to_blen_seg_pairs[isostate].append(blen_seg_pair)
    seg_to_text_label = {}
    for isostate, blen_seg_pairs in isostate_to_blen_seg_pairs.items():
        blen, seg = max(blen_seg_pairs)
        state = seg_to_state[seg]
        seg_to_text_label[seg] = str(state)

    # draw the image
    surface = cairo.SVGSurface('out.svg', screenwidth, screenheight)
    ctx = cairo.Context(surface)
    ctx.translate(screenwidth/2, screenheight/2)
    for va, vb in G.edges():
        seg = G[va][vb]['segment']
        rgb = parity_to_rgb[G[va][vb]['parity']]
        pa = (vertex_to_2d[va] - mid) * scale
        pb = (vertex_to_2d[vb] - mid) * scale
        ctx.set_source_rgb(*rgb)
        ctx.move_to(pa[0], pa[1])
        ctx.line_to(pb[0], pb[1])
        ctx.stroke()
        text_label = seg_to_text_label.get(seg, None)
        if text_label is not None:
            pmid = (pa + pb) / 2.0
            text_extents = ctx.text_extents(text_label)
            x_bearing, y_bearing, width, height = text_extents[:4]
            ctx.set_source_rgb(*black_rgb)
            ctx.move_to(
                    pmid[0] - width / 2.0 - x_bearing,
                    pmid[1] - height / 2.0 - y_bearing,
                    )
            ctx.show_text(text_label)

    # finish drawing the image and create the image file
    ctx = None
    surface.finish()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--history', default='history.db',
            help='sampled state history as an sqlite3 database file')
    parser.add_argument('--layout', default='layout.db',
            help='2d node layout as an sqlite3 database file')
    parser.add_argument('--segmentation', default='segmentation.db',
            help='isostate segmentation as an sqlite3 database file')
    #parser.add_argument('--imageformat', choices=('png', 'svg', 'pdf'),
            #help='choose among image formats')
    parser.add_argument('-o', '--outfile', default='out.png',
            help='create this image file')
    main(parser.parse_args())

