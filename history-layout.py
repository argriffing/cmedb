"""
Associate 2d coordinates with each vertex of a sampled history.

This script should possibly be redesigned
for more suitable input and output.
"""

import argparse
import sqlite3

import numpy as np
import networkx as nx


def main(args):

    # open the history db file
    # read the association of node index pairs to edge indices
    # read the association of edge indices to branch lengths
    conn = sqlite3.connect(args.infile)
    cursor = conn.cursor()
    cursor.execute('select va, vb, blen from history')
    print list(cursor)
    conn.close()

    # create or open the database
    #conn = sqlite3.connect('rate.matrix.db')
    #cursor = conn.cursor()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', required=True,
            help='sampled state history as an sqlite3 database file')
    parser.add_argument('-o', '--outfile', required=True,
            help='create this sqlite3 database file')
    main(parser.parse_args())

