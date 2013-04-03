"""
For the blinking process, partition codon states by amino acid.

The input is a genetic code, and the output is a state partition.
"""

import argparse
import sqlite3


def main(args):

    # read the genetic code
    conn = sqlite3.connect(args.code)
    cursor = conn.cursor()
    cursor.execute(
        "select state, residue, codon from code "
        "where residue <> 'Stop' "
        "order by state")
    genetic_code = list(cursor)
    conn.close()

    # get a list of distinct residues
    states, residues, codons = zip(*genetic_code)
    residues = sorted(set(residues))
    residue_to_part = dict((r, i) for i, r in enumerate(residues))

    # open the database file for the patterns
    conn = sqlite3.connect(args.outfile)
    cursor = conn.cursor()

    # initialize the pattern multiplicities table
    s = (
            'create table if not exists partition ('
            'state integer, '
            'part integer, '
            'primary key (state))')
    cursor.execute(s)
    conn.commit()

    # write the partition information into the table
    for state, residue, codon in genetic_code:
        part = residue_to_part[residue]
        s = 'insert into partition values (?, ?)'
        t = (state, part)
        cursor.execute(s, t)
    conn.commit()

    # close the database file for the patterns
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--code', default='universal.code.db',
            help='input genetic code in sqlite3 format')
    parser.add_argument('-o', '--outfile', default='codon.partition.db',
            help='create this sqlite3 database file')
    main(parser.parse_args())

