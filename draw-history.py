"""
Create an image file visualization of a state history on a tree using pycairo.

The layout should be equal angle without much complicatedness.
For example, do not bother optimizing the rotation yet.
On the other hand, take care to draw branch lengths proportionally.
More complicated layouts like equal daylight or layouts
that try to iteratively maximimize a layout goodness are not implemented.
"""

#FIXME add some code...
#FIXME break into two scripts.
#FIXME the first script (layout) will compute the vertex coordinates
#FIXME the second script (rendering) will draw the tree to an image file

def main(args):
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tree', required=True,
            help='unrooted tree as an sqlite3 database file')
    parser.add_argument('--rates', required=True,
            help=('continuous-time Markov substitution process '
                'as an sqlite3 database file'))
    parser.add_argument('--root', type=int,
            help='root node index')
    parser.add_argument('--nsamples', type=pos_int, default=5,
            help='sample this many histories')
    parser.add_argument('-o', '--outfile', required=True,
            help='create this sqlite3 database file')
    parser.add_argument('--table', required=True,
            help='name of table to create in the new database')
    main(parser.parse_args())

