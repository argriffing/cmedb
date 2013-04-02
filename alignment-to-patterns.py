"""
Extract site patterns from an alignment.
"""

import argparse
import sqlite3
from collections import defaultdict


def main(args):

    # read the alignment from the database
    conn = sqlite3.connect(args.alignment)
    cursor = conn.cursor()
    cursor.execute(
            'select offset, taxon, state from alignment '
            'order by offset, taxon')
    offset_taxon_state = list(cursor)
    conn.close()

    # reduce the assignments to patterns with multiplicities
    offset_to_pairs = defaultdict(list)
    for offset, taxon, state in offset_taxon_state:
        offset_to_pairs[offset].append((taxon, state))
    pairs_to_count = defaultdict(int)
    for pairs in offset_to_pairs.values():
        pairs_to_count[tuple(pairs)] += 1
    pairs_count = [x for x in pairs_to_count.items()]

    # open the database file for the patterns
    conn = sqlite3.connect(args.outfile)
    cursor = conn.cursor()

    # initialize the pattern multiplicities table
    s = (
            'create table if not exists multiplicities ('
            'pattern integer, '
            'multiplicity integer, '
            'primary key (pattern))')
    cursor.execute(s)
    conn.commit()

    # initialize the pattern definition table
    s = (
            'create table if not exists patterns ('
            'pattern integer, '
            'taxon integer, '
            'state integer, '
            'primary key (pattern, taxon))')
    cursor.execute(s)
    conn.commit()

    # write the pattern multiplicities into a table
    for pattern_index, (pairs, count) in enumerate(pairs_count):
        s = 'insert into multiplicities values (?, ?)'
        t = (pattern_index, count)
        cursor.execute(s, t)
    conn.commit()

    # write the pattern definitions into a table
    for pattern_index, (pairs, count) in enumerate(pairs_count):
        for taxon, state in pairs:
            s = 'insert into patterns values (?, ?, ?)'
            t = (pattern_index, taxon, state)
            cursor.execute(s, t)
    conn.commit()

    # close the database file for the patterns
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--alignment', default='alignment.db',
            help='input i.i.d. aligned states as an sqlite3 database file')
    parser.add_argument('-o', '--outfile', default='patterns.db',
            help='output alignment patterns as an sqlite3 database file')
    main(parser.parse_args())

