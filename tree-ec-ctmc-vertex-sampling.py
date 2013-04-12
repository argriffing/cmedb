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
The ec in the script name refers to "endpoint conditioned"
which means "leaf conditioned" in the case of a tree.
More generally it is "boundary conditioned."
The ctmc in the script name refers to "continuous time Markov chain."
This script in particular assumes that the Markov chain is time-reversible.
"""

# An experimental way to organize scripts,
# possibly to be used with galaxy.
g_tags = {'tree', 'sampling', 'ctmc', 'reversible'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tree', default='brown.tree.db',
            help='unrooted tree as an sqlite3 database file')
    parser.add_argument('--rates', default='rate.matrix.db',
            help=('continuous-time Markov substitution process '
                'as an sqlite3 database file'))
    parser.add_argument('--site', default='site.db',
            help='aligned states at the leaves of the tree')
    parser.add_argument('--outfile', default='augmented.alignment.db',
            help='augmented alignment at all nodes of the tree')
    main(parser.parse_args())

