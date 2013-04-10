"""
We are making the blinking process complicated by adding heterogeneous rates.
"""

import random
import argparse
import sqlite3

import cmedbutil


def main(args):

    # Read the partition from a database file.
    # The hidden blinking process controls the state transition
    # of the primary process according to this partition
    # of primary process states.
    conn = sqlite3.connect(args.partition)
    cursor = conn.cursor()
    cursor.execute('select part from partition')
    parts = set(t[0] for t in cursor)
    conn.close()

    # open the database file for the blink rates
    conn = sqlite3.connect(args.outfile)
    cursor = conn.cursor()

    # initialize the blink rates table
    s = (
            'create table if not exists rates ('
            'part integer, '
            'rate_on real, '
            'rate_off real, '
            'primary key (part))')
    cursor.execute(s)
    conn.commit()

    # write the partition information into the table
    mu = 1.0
    rate_on = random.expovariate(1 / mu)
    rate_off = random.expovariate(1 / mu)
    for part in parts:
        s = 'insert into rates values (?, ?, ?)'
        
        # define the rate on
        if args.rate_on is not None:
            rate_on = args.rate_on
        elif not args.shared_rates:
            rate_on = random.expovariate(1 / mu)

        # define the rate off
        if args.rate_off is not None:
            rate_off = args.rate_off
        elif not args.shared_rates:
            rate_off = random.expovariate(1 / mu)

        # add a row to the database
        t = (part, rate_on, rate_off)
        cursor.execute(s, t)
    conn.commit()

    # close the database file for the patterns
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--partition', default='partition.db',
            help=('a partition of the primary states '
                'as an sqlite3 database file'))
    parser.add_argument('--outfile', default='blink.rates.db',
            help='create this sqlite3 database file')
    parser.add_argument('--shared-rates', action='store_true',
            help='all blink states share the same rates')
    parser.add_argument('--rate-on', type=cmedbutil.nonneg_float,
            help='use a fixed blink rate on')
    parser.add_argument('--rate-off', type=cmedbutil.nonneg_float,
            help='use a fixed blink rate off')
    main(parser.parse_args())

