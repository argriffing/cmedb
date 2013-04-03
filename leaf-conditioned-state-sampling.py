"""
Given states at the leaves of an unrooted tree, sample internal vertex states.

This sampling requires a rate matrix
for which the expm is computed on each branch of the tree.
A Felsenstein-like dynamic programming algorithm is used.
Note that the input leaf alignment is assumed to be
separated into distinct column patterns and their multiplicities,
whereas the output alignment
does not pay attention pattern distinctness or multiplicity.
Also the leaf state alignment is unnecessarily
duplicated in the output alignment.
"""

#XXX obviously not finished

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tree', default='brown.tree.db',
            help='unrooted tree as an sqlite3 database file')
    parser.add_argument('--rates', default='rate.matrix.db',
            help=('continuous-time Markov substitution process '
                'as an sqlite3 database file'))
    parser.add_argument('--leaf-patterns', default='patterns.db',
            help='aligned states at the leaves of the tree')
    parser.add_argument('--outfile', default='augmented.alignment.db',
            help='augmented alignment at all nodes of the tree')
    main(parser.parse_args())

