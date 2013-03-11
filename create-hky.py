"""
This script creates a nucleotide rate matrix sqlite database.
"""

import argparse
import sqlite3

import numpy as np

def main(args):

    # check the command line arguments
    distn = np.array([args.A, args.C, args.G, args.T])
    if not np.allclose(np.sum(distn), 1):
        raise Exception('nucleotide distribution should sum to 1.0')
    if np.any(distn < 0):
        raise Exception('nucleotide probabilities should be non-negative')

    # create or open the database
    conn = sqlite3.connect('rate.matrix.db')
    cursor = conn.cursor()

    # create the tables
    cursor.execute(
            'create table rates ('
            'source integer, sink integer, rate real, '
            'primary key (source, sink))')
    cursor.execute(
            'create table states ('
            'state integer, name text, '
            'primary key (state))')
    cursor.execute(
            'create table distn ('
            'state integer, prob real, '
            'primary key (state))')
    conn.commit()

    # define the nucleotide states
    states = 'ACGT'

    # populate the rate matrix table
    for i, si in enumerate(states):
        for j, sj in enumerate(states):
            if i == j:
                rate = 0
            elif si+sj in ('AG', 'GA', 'CT', 'TC'):
                rate = args.kappa * distn[j]
            else:
                rate = distn[j]
            if rate > 0:
                cursor.execute(
                        'insert into rates values (?, ?, ?)',
                        (i, j, rate),
                        )
    conn.commit()

    # populate the states table
    for pair in enumerate(states):
        cursor.execute('insert into states values (?, ?)', pair)
    conn.commit()

    # populate the state distribution table
    for pair in enumerate(distn):
        cursor.execute('insert into distn values (?, ?)', pair)
    conn.commit()

    # close the database connection
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-A', default=0.25, type=float,
            help='probability of nucleotide A')
    parser.add_argument('-C', default=0.25, type=float,
            help='probability of nucleotide C')
    parser.add_argument('-G', default=0.25, type=float,
            help='probability of nucleotide G')
    parser.add_argument('-T', default=0.25, type=float,
            help='probability of nucleotide T')
    parser.add_argument('-k', '--kappa', default=2, type=float,
            help='transition/transversion ratio')
    main(parser.parse_args())

