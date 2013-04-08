"""
Compute Monte Carlo likelihood ratio given a couple of log likeihood files.
"""

import argparse
import math
import sqlite3

def main(args):

    # define the sql command to read the log likelihoods
    s = 'select history, log_likelihood from log_likelihoods'

    # read numerator log likelihoods
    conn = sqlite3.connect(args.numerator_log_likelihoods)
    cursor = conn.cursor()
    history_to_numerator = dict(cursor.execute(s))
    conn.close()

    # read numerator log likelihoods
    conn = sqlite3.connect(args.denominator_log_likelihoods)
    cursor = conn.cursor()
    history_to_denominator = dict(cursor.execute(s))
    conn.close()

    # verify that the history indices are the same
    if set(history_to_numerator) != set(history_to_denominator):
        raise Exception

    # get the expected likelihood ratio
    history_indices = set(history_to_numerator)
    nhistories = len(history_indices)
    accum = 0.0
    for history in history_indices:
        log_num = history_to_numerator[history]
        log_den = history_to_denominator[history]
        accum += math.exp(log_num - log_den)
    mclr = accum / nhistories

    # report the expected likelihood ratio
    print 'Monte Carlo likelihood ratio:', mclr


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--numerator-log-likelihoods',
            default='log.likelihoods.db',
            help='numerator log likelihoods as an sqlite3 database file')
    parser.add_argument('--denominator-log-likelihoods',
            default='log.likelihoods.db',
            help='denominator log likelihoods as an sqlite3 database file')
    main(parser.parse_args())
